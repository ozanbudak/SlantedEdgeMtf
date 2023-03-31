using MathNet.Numerics;
using Numpy;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Windows;

namespace SlantedEdgeMtf
{
    public static class MtfAlgorithm
    {
        public static (double[] esf, double[] psf, double nfreq, double[] freqs, double[] mtf, double slope, double eff) CalculateMtf(Bitmap bitmap, double sampling_interval_mm, int binning_factor)
        {
            var nlin = bitmap.Height;

            var npix = bitmap.Width;

            double del = sampling_interval_mm;

            int nbin = binning_factor;

            var a = BitmapToArray(bitmap);

            var test = ContrastTest(a);

            if (test < 0.2)
            {
                MessageBox.Show("Edge contrast is less that 20%, this can lead to high error in the SFR measurement.", "WARNING", MessageBoxButton.OK, MessageBoxImage.Warning);
            }

            var fil1 = new double[2] { -0.5, 0.5 };

            var fil2 = new double[3] { -0.5, 0, 0.5 };

            var win1 = CreateHammingWindow(npix, (npix + 1) / 2);

            var c = Derivative(a, nlin, npix, fil1);

            var loc = new double[nlin];

            for (int y = 0; y < nlin; y++)
            {
                var tempArray = new double[npix];

                for (int x = 0; x < npix; x++)
                {
                    tempArray[x] = c[y, x] * win1[x];
                }

                loc[y] = Math.Round(Centroid(tempArray) - 0.5, 4);
            }

            var fitme = FindEdge(loc, nlin);

            var place = new double[nlin];

            for (int y = 0; y < nlin; y++)
            {
                place[y] = fitme.intercept + fitme.slope * y;
                var win2 = CreateHammingWindow(npix, place[y]);

                var tempArray = new double[npix];

                for (int x = 0; x < npix; x++)
                {
                    tempArray[x] = c[y, x] * win2[x];
                }

                loc[y] = Math.Round(Centroid(tempArray) - 0.5, 4);
            }

            fitme = FindEdge(loc, nlin);

            var vslope = fitme.slope;

            var slope_deg = Math.Round(180 * Math.Round(Math.Atan(Math.Abs(vslope)), 4) / Math.Round(Math.PI, 4), 4);

            if (slope_deg < 3.5)
            {
                MessageBox.Show($" High slope, {slope_deg:N2} degrees.", "WARNING", MessageBoxButton.OK, MessageBoxImage.Warning);
            }

            var delimage = del;

            var delfac = Math.Round(Math.Cos(Math.Atan(vslope)), 4);

            del = Math.Round(del * delfac, 4);

            var del2 = del / nbin;

            var nn = (int)Math.Floor((double)npix * nbin);

            int nn2 = (int)Math.Floor(nn / 2.0) + 1;

            var dcorr = Fir2Fix(nn2, 3);

            var freqs = new double[nn];

            for (int i = 0; i < nn; i++)
            {
                freqs[i] = nbin * i / (del * nn);
            }

            var freqlim = 1;

            if (nbin == 1)
                freqlim = 2;

            var nn2out = Math.Round((double)nn2 * freqlim / 2, 0);

            var nfreq = nn / (2 * delimage * nn);

            var win = CreateHammingWindow(nbin * npix, (double)(nbin * npix + 1) / 2);

            var esf = new double[nn];

            var point = Project(a, fitme.slope, nbin);

            esf = point;

            var point2d = Make2DArray(point, 1, point.Length);

            var c1 = Derivative(point, fil2);

            var psf = c1;

            var mid = Centroid(c1);

            var temp = Centroid(c1, (int)Math.Round(mid));

            c1 = temp;

            for (int x = 0; x < c1.Length; x++)
            {
                c1[x] = c1[x] * win[x];
            }

            var psf_ = c1;

            var fftc1 = np.abs(np.fft.fft_(c1)).GetData<double>();

            var mtf = fftc1.Take(nn / 2).Select(x => x / fftc1[0]).ToArray();

            for (int i = 0; i < mtf.Length; i++)
            {
                mtf[i] = mtf[i] * dcorr[i];
            }

            var dat = new double[(int)nn2out, 2];

            for (int i = 0; i < nn2out; i++)
            {
                dat[i, 0] = freqs[i];
                dat[i, 1] = mtf[i];
            }

            var val = new double[2] { 0.1, 0.5 };

            var (eff, freqval, sfrval) = SampleEff(dat, val, delimage);

            return (esf, psf_, nfreq, freqs.Take((int)nn2out).ToArray(), mtf.Take((int)nn2out).ToArray(), slope_deg, eff);
        }

        public static double[,] BitmapToArray(Bitmap bitmap)
        {
            var bitmap_array = new double[bitmap.Height, bitmap.Width];

            for (int x = 0; x < bitmap.Width; x++)
            {
                for (int y = 0; y < bitmap.Height; y++)
                {
                    bitmap_array[y, x] = bitmap.GetPixel(x, y).R;
                }
            }

            return bitmap_array;
        }

        public static double ContrastTest(double[,] bitmap_array)
        {
            double tleft = 0;

            for (int y = 0; y < bitmap_array.GetLength(0); y++)
            {
                for (int x = 0; x < 5; x++)
                {
                    tleft += bitmap_array[y, x];
                }
            }

            double tright = 0;

            for (int y = 0; y < bitmap_array.GetLength(0); y++)
            {
                for (int x = bitmap_array.GetLength(1) - 6; x < bitmap_array.GetLength(1); x++)
                {
                    tright += bitmap_array[y, x];
                }
            }

            return Math.Abs((tleft - tright) / (tleft + tright));
        }

        public static double[] CreateHammingWindow(int n, double mid)
        {
            double[] data = new double[n];

            var wid1 = mid - 1;

            var wid2 = n - mid;

            var wid = Math.Max(wid1, wid2);

            double pie = Math.Round(Math.PI, 4);

            for (int i = 0; i < n; i++)
            {
                double arg = (i + 1) - mid;
                data[i] = Math.Cos(pie * arg / wid);
            }

            for (int i = 0; i < n; i++)
            {
                data[i] = Math.Round(0.54 + (0.46 * data[i]), 4);
            }

            return data;
        }

        public static double[,] Derivative(double[,] a, int nlin, int npix, double[] fil)
        {
            double[,] b = new double[nlin, npix];

            int nn = fil.Length;

            for (int i = 0; i < nlin; i++)
            {
                var array = new double[npix];

                for (int k = 0; k < npix; k++)
                {
                    array[k] = a[i, k];
                }

                var conv = Convolve(array, fil, true).Select(x => Math.Round(x, 4)).ToArray();

                for (int conv_index = 0; conv_index < conv.Length - 1; conv_index++)
                {
                    conv[conv_index] = conv[conv_index + 1];
                }

                conv[0] = conv[1];

                conv[^1] = conv[^2];


                for (int j = 0; j < conv.Length; j++)
                {
                    b[i, j] = conv[j];
                }
            }

            return b;
        }

        public static double[] Derivative(double[] a, double[] fil)
        {
            var conv = Convolve(a, fil, true).ToArray();

            conv[0] = conv[1];

            conv[^1] = conv[^2];

            return conv;
        }

        public static double Centroid(double[] x)
        {
            double n = x.Length;

            double sumx = x.Sum();

            double loc;

            if (sumx < 1e-4)
            {
                loc = 0;
            }
            else
            {
                double sumNx = 0;
                for (int i = 0; i < n; i++)
                {
                    sumNx += (i + 1) * x[i];
                }
                loc = sumNx / sumx;
            }

            return loc;
        }

        public static (double slope, double intercept) FindEdge(double[] cent, int nlin)
        {
            var index = Enumerable.Range(0, nlin).Select(x => (double)x).ToArray();

            var (intercept, slope) = Fit.Line(index, cent);

            return (slope, intercept);
        }

        public static double[] Fir2Fix(int n, int m)
        {
            double[] correct = new double[n];

            m--;

            double scale = 1.0;

            for (int i = 1; i < n; i++)
            {
                correct[i] = Math.Abs((Math.PI * i * m / (2.0 * (n + 1.0))) / Math.Sin(Math.PI * i * m / (2.0 * (n + 1.0))));
                correct[i] = Math.Round(1.0 + scale * (correct[i] - 1.0), 4);

                if (correct[i] > 10.0)
                {
                    correct[i] = 10.0;
                }
            }

            correct[0] = 1;

            return correct;
        }

        public static double[] Project(double[,] bb, double slope, double fac)
        {
            var nlin = bb.GetLength(0);

            var npix = bb.GetLength(1);

            var nn = npix * fac;

            slope = 1 / slope;

            var offset = Math.Round(fac * (0 - (nlin - 1) / slope));

            var del = Math.Abs(offset);

            if (offset > 0)
                offset = 0;

            double[,] barray = new double[2, (int)(nn + del + 100)];

            for (int n = 1; n <= npix; n++)
            {
                for (int m = 1; m <= nlin; m++)
                {
                    double x = n - 1;

                    double y = m - 1;

                    int ling = (int)(Math.Ceiling((x - y / slope) * fac) + 1 - offset);

                    barray[0, ling] += 1;

                    barray[1, ling] += bb[m - 1, n - 1];
                }
            }

            var point = new double[(int)nn];

            var start = 1 + Math.Round(0.5 * del);

            int nz = 0;

            for (int i = (int)start; i < start + nn; i++)
            {
                if (barray[0, i] == 0)
                {
                    nz++;

                    MessageBox.Show("Zero count(s) found during projection binning. The edge angle may be large, or you may need more lines of data. Execution will continue.", "WARNING", MessageBoxButton.OK, MessageBoxImage.Warning);

                    if (i == 1)
                    {
                        barray[0, i] = barray[0, i + 1];
                    }
                    else
                    {
                        barray[0, i] = (barray[0, i - 1] + barray[0, i + 1]) / 2;
                    }
                }
            }

            for (int i = 0; i < nn; i++)
            {
                point[i] = barray[1, i + (int)start] / barray[0, i + (int)start];
            }

            return point;
        }

        public static double[] Centroid(double[] a, int center)
        {
            int n = a.Length;

            double[] b = new double[n];

            int mid = (int)(Math.Round((double)(n + 1) / 2, 0, MidpointRounding.ToPositiveInfinity));

            int del = center - mid;

            if (del > 0)
            {
                for (int i = 0; i < n - del; i++)
                {
                    b[i] = a[i + del];
                }
            }
            else if (del < 0)
            {
                for (int i = -del; i < n; i++)
                {
                    b[i] = a[i + del];
                }
            }
            else
            {
                b = a;
            }

            return b;
        }

        public static (double eff, double[] freqval, double[] sfrval) SampleEff(double[,] dat, double[] val, double del)
        {
            var mmin = val.Min();

            var mindex = val.ToList().IndexOf(mmin);

            if (mmin > 0.1)
            {
                MessageBox.Show($"WARNING: sampling efficiency is based on SFR {mmin}", "WARNING", MessageBoxButton.OK, MessageBoxImage.Warning);
            }

            var delf = dat[1, 0] + 1e-6;

            var hs = 0.5 / del;

            List<int> xList = new();

            for (int i = 0; i < dat.GetLength(0); i++)
            {
                if (dat[i, 0] > 1.1 * hs)
                {
                    xList.Add(i);
                }
            }

            int[] x = xList.ToArray();

            double imax;

            double imaxx;

            if (x.Length == 0)
            {
                imax = dat.GetLength(0);
                imaxx = imax;
            }
            else
            {
                List<int> xxList = new();

                for (int i = 0; i < dat.GetLength(0); i++)
                {
                    if (dat[i, 0] > hs - delf)
                    {
                        xxList.Add(i);
                    }
                }
                int[] xx = xxList.ToArray();

                imax = x[0];

                imaxx = xx[0];

                double[,] newDat = new double[(int)imax + 1, dat.GetLength(1)];

                for (int i = 0; i < imax + 1; i++)
                {
                    for (int j = 0; j < dat.GetLength(1); j++)
                    {
                        newDat[i, j] = dat[i, j];
                    }
                }

                dat = newDat;
            }

            var nval = val.Length;

            double[] freqval = new double[nval];

            double[] sfrval = new double[nval];

            for (int v = 0; v < nval; v++)
            {
                var sftfreq = FindFrequency(dat, val[v], imax);

                freqval[v] = sftfreq.freqval[0];

                sfrval[v] = sftfreq.sfrval[0];
            }

            double eff = Math.Min(Math.Round(100 * freqval[mindex] / dat[(int)imaxx, 0]), 100);

            return (eff, freqval, sfrval);
        }

        public static (double[] freqval, double[] sfrval) FindFrequency(double[,] dat, double val, double imax)
        {
            int n = dat.GetLength(0);

            int m = dat.GetLength(1);

            int nc = m - 1;

            double[] freqval = new double[nc];

            double[] sfrval = new double[nc];

            var maxf = dat[(int)imax, 0];

            var test = new double[n];

            for (int c = 0; c < nc; c++)
            {
                for (int i = 0; i < n; i++)
                {
                    test[i] = dat[i, c + 1] - val;
                }

                List<int> xList = new();

                for (int i = 0; i < test.Length; i++)
                {
                    if (test[i] < 0)
                    {
                        xList.Add(i);
                    }
                }
                int[] x = xList.ToArray();

                var x1 = x[0];

                var sval = dat[x1 - 1, c + 1];

                var s = dat[x1 - 1, 0];

                var y = dat[x1 - 1, c + 1];

                var y2 = dat[x1, c + 1];

                var slope = (y2 - y) / dat[1, 0];

                var dely = test[x1 - 1];

                s -= dely / slope;

                sval -= dely;

                if (s > maxf)
                {
                    s = maxf;
                    sval = dat[(int)imax, c + 1];
                }

                freqval[c] = s;

                sfrval[c] = sval;
            }

            return (freqval, sfrval);
        }

        private static T[,] Make2DArray<T>(T[] input, int height, int width)
        {
            T[,] output = new T[height, width];

            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    output[i, j] = input[i * width + j];
                }
            }

            return output;
        }

        public static double[] Convolve(this double[] a, double[] kernel, bool trim)
        {
            double[] result;

            int m = (int)Math.Ceiling(kernel.Length / 2.0);

            if (trim)
            {
                result = new double[a.Length];

                for (int i = 0; i < result.Length; i++)
                {
                    result[i] = 0;
                    for (int j = 0; j < kernel.Length; j++)
                    {
                        int k = i - j + m - 1;
                        if (k >= 0 && k < a.Length)
                            result[i] += a[k] * kernel[j];
                    }
                }
            }
            else
            {
                result = new double[a.Length + m];

                for (int i = 0; i < result.Length; i++)
                {
                    result[i] = 0;
                    for (int j = 0; j < kernel.Length; j++)
                    {
                        int k = i - j;
                        if (k >= 0 && k < a.Length)
                            result[i] += a[k] * kernel[j];
                    }
                }
            }

            return result;
        }
    }
}
