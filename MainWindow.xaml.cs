using MathNet.Numerics.Interpolation;
using Microsoft.Win32;
using System;
using System.Drawing.Imaging;
using System.Drawing;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media.Imaging;
using Rectangle = System.Drawing.Rectangle;
using Color = System.Drawing.Color;

namespace SlantedEdgeMtf
{

    public partial class MainWindow : Window
    {
        public Rectangle RoiRect { get; set; } = new(0, 0, 256, 600);
        public double MicronPerPixel { get; set; } = 2;
        public double FirstReadOutLpmm { get; set; } = 2.5;
        public double SecondReadOutLpmm { get; set; } = 7.5;
        public double ThirdReadOutLpmm { get; set; } = 15;
        public double FourthReadOutLpmm { get; set; } = 30;
        public Bitmap? ImageBitmap { get; set; }
        public Bitmap? ImageBitmapRoi { get; set; }

        public MainWindow()
        {
            InitializeComponent();
        }

        public void ProcessImage(Bitmap bitmap)
        {
            try
            {
                if (bitmap != null)
                {
                    var (esf, psf, nfreq, freqs, mtf, slope, eff) = MtfAlgorithm.CalculateMtf(bitmap, MicronPerPixel / 1000, 4);

                    EsfPlot.Visibility = Visibility.Visible;

                    EsfPlot.Plot.Clear();

                    EsfPlot.Plot.Title("ESF", fontName: "Arial");

                    EsfPlot.Plot.AddSignal(esf.Reverse().ToArray(), 1, Color.Black, "Edge Spread Function");

                    EsfPlot.Plot.AddText($"{bitmap.Width} x {bitmap.Height} pixels (WxH)", .1 * nfreq, 180, size: 12, color: Color.Black);

                    EsfPlot.Plot.AddText($"{MicronPerPixel} um per pixel", .1 * nfreq, 160, size: 12, color: Color.Black);

                    EsfPlot.Plot.AddText($"Edge {slope:N2} °", .1 * nfreq, 140, size: 12, color: Color.Black);

                    EsfPlot.Plot.XLabel("Pixels (4x oversampled)");

                    EsfPlot.Plot.YLabel("Edge Pixel Profile");

                    EsfPlot.Refresh();

                    MtfPlot.Plot.Clear();

                    MtfPlot.Plot.Title("MTF", fontName: "Arial");

                    MtfPlot.Plot.AddSignalXY(freqs, mtf, Color.Black, "Spatial Frequency Response");

                    MtfPlot.Plot.AddAnnotation($"Sampling Efficiency {eff}%", -10, 10);

                    MtfPlot.Plot.AddAnnotation($"Slope {slope:N2} °", -10, 40);

                    MtfPlot.Plot.AddText("Half-sampling", .9 * nfreq, .125, size: 13, color: Color.Black);

                    IInterpolation interpolation = LinearSpline.Interpolate(freqs, mtf);

                    double interpolatedValue = interpolation.Interpolate(FirstReadOutLpmm);

                    double interpolatedValue2 = interpolation.Interpolate(SecondReadOutLpmm);

                    double interpolatedValue3 = interpolation.Interpolate(ThirdReadOutLpmm);

                    double interpolatedValue4 = interpolation.Interpolate(FourthReadOutLpmm);

                    MtfPlot.Plot.AddText($"MTF @ {FirstReadOutLpmm}lp/mm = {interpolatedValue:N3}", .6 * nfreq, 1, size: 12, color: Color.Black);

                    MtfPlot.Plot.AddText($"MTF @ {SecondReadOutLpmm}lp/mm = {interpolatedValue2:N3}", .6 * nfreq, 0.9, size: 12, color: Color.Black);

                    MtfPlot.Plot.AddText($"MTF @ {ThirdReadOutLpmm}lp/mm = {interpolatedValue3:N3}", .6 * nfreq, 0.8, size: 12, color: Color.Black);

                    MtfPlot.Plot.AddText($"MTF @ {FourthReadOutLpmm}lp/mm = {interpolatedValue4:N3}", .6 * nfreq, 0.7, size: 12, color: Color.Black);

                    MtfPlot.Plot.XLabel($"Frequency lp/mm, ({MicronPerPixel} um/Pixel)");

                    MtfPlot.Plot.YLabel("SFR");

                    MtfPlot.Plot.SetAxisLimitsY(0, 1.2);

                    var vLine1 = MtfPlot.Plot.AddVerticalLine(nfreq, width: 1);
                    vLine1.Max = 0.05;
                    vLine1.Color = Color.Blue;

                    MtfPlot.Refresh();

                    PsfPlot.Plot.Clear();

                    PsfPlot.Plot.Title("PSF", fontName: "Arial");

                    PsfPlot.Plot.AddSignal(psf, 1, Color.Black, "Point Spread Function");

                    PsfPlot.Refresh();
                }
            }
            catch (Exception exception)
            {
                MessageBox.Show(exception.Message, Math.Abs(exception.Message.GetHashCode()).ToString(CultureInfo.CurrentCulture));
            }
        }

        private void RoiAppyButton_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                var roiX = Convert.ToInt32(RoiXPositionTextBox.Text);

                var roiY = Convert.ToInt32(RoiYPositionTextBox.Text);

                var roiWidth = Convert.ToInt32(RoiWidthTextBox.Text);

                var roiHeight = Convert.ToInt32(RoiHeightTextBox.Text);

                if (ImageBitmap is not null)
                {
                    if (roiX > ImageBitmap.Width || roiY > ImageBitmap.Height)
                    {
                        MessageBox.Show("ROI width and heigth cannot be higher than image width and height.");
                        return;
                    }

                    if (roiWidth > ImageBitmap.Width - roiX)
                    {
                        MessageBox.Show($"ROI width cannot be larger than total width of image is starting from {roiX}. Max Width is {ImageBitmap.Width - roiX}");
                        return;
                    }

                    if (roiHeight > ImageBitmap.Height - roiY)
                    {
                        MessageBox.Show($"ROI height cannot be larger than total height of image is starting from {roiY}. Max Width is {ImageBitmap.Height - roiY}");
                        return;
                    }

                    RoiRect = new(roiX, roiY, roiWidth, roiHeight);

                    var newBitmap = ImageBitmap.Clone(new Rectangle(roiX, roiY, roiWidth, roiHeight), ImageBitmap.PixelFormat);

                    ImageBitmapRoi = newBitmap;

                    RoiImage.Source = BitmapToImageSource(ImageBitmapRoi);

                    ProcessImage(ImageBitmapRoi);
                }
            }
            catch (Exception exception)
            {
                MessageBox.Show(exception.Message, Math.Abs(exception.Message.GetHashCode()).ToString(CultureInfo.CurrentCulture));
            }
        }

        private void DisplayApplyButton_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                MicronPerPixel = Convert.ToDouble(MicronsPerPixelTextBox.Text);

                FirstReadOutLpmm = Convert.ToDouble(FirstReadOutTextBox.Text);

                SecondReadOutLpmm = Convert.ToDouble(SecondReadOutTextBox.Text);

                ThirdReadOutLpmm = Convert.ToDouble(ThirdReadOutTextBox.Text);

                FourthReadOutLpmm = Convert.ToDouble(FourthReadOutTextBox.Text);

                if (ImageBitmap != null)
                {
                    ProcessImage(ImageBitmap);
                }
            }
            catch (Exception exception)
            {
                MessageBox.Show(exception.Message, Math.Abs(exception.Message.GetHashCode()).ToString(CultureInfo.CurrentCulture));
            }
        }

        private void SelectImageButton_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                var openFileDialog = new OpenFileDialog
                {
                    InitialDirectory = @"C:\",
                    Title = "Browse PNG Files",

                    CheckFileExists = true,
                    CheckPathExists = true,

                    DefaultExt = "png",
                    Filter = "png files (*.png)|*.png",
                    FilterIndex = 2,
                    RestoreDirectory = true,

                    ReadOnlyChecked = true,
                    ShowReadOnly = true
                };

                if (openFileDialog.ShowDialog() == true)
                {
                    var fileName = openFileDialog.FileName;

                    var bmp = new Bitmap(fileName);

                    ImageBitmap = MakeGrayscale3(bmp);

                    RoiImage.Source = BitmapToImageSource(ImageBitmap);

                    RoiRect = new Rectangle(0, 0, ImageBitmap.Width, ImageBitmap.Height);

                    RoiWidthTextBox.Text = ImageBitmap.Width.ToString();

                    RoiHeightTextBox.Text = ImageBitmap.Height.ToString();

                    ProcessImage(ImageBitmap);
                }
            }
            catch (Exception ex)
            {
                MessageBox.Show($"WARNING: Unknown error, {ex.Message}.", "WARNING", MessageBoxButton.OK, MessageBoxImage.Warning);
            }
        }

        private void DecimalValidation_PreviewTextInput(object sender, TextCompositionEventArgs e)
        {
            Regex regex = new("^[-][,][0-9]+$|^[-][0-9]*[,]{0,1}[0-9]*$|^[0-9]*[,]{0,1}[0-9]*$");
            e.Handled = !regex.IsMatch(((TextBox)sender).Text.Insert(((TextBox)sender).SelectionStart, e.Text));
        }

        private void IntegerValidation_PreviewTextInput(object sender, TextCompositionEventArgs e)
        {
            var regex = new Regex("[^0-9]+");
            e.Handled = regex.IsMatch(e.Text);
        }

        private void CloseButton_Click(object sender, RoutedEventArgs e)
        {
            Close();
        }

        private void HeaderGrid_MouseDown(object sender, MouseButtonEventArgs e)
        {
            if (e.ChangedButton == MouseButton.Left)
                DragMove();
        }

        public static BitmapImage BitmapToImageSource(Bitmap bitmap)
        {
            using MemoryStream memory = new();

            bitmap.Save(memory, ImageFormat.Bmp);

            memory.Position = 0;

            BitmapImage bitmapImage = new();

            bitmapImage.BeginInit();

            bitmapImage.StreamSource = memory;

            bitmapImage.CacheOption = BitmapCacheOption.OnLoad;

            bitmapImage.EndInit();

            return bitmapImage;
        }

        public static Bitmap MakeGrayscale3(Bitmap original)
        {

            Bitmap newBitmap = new(original.Width, original.Height);

            using (Graphics g = Graphics.FromImage(newBitmap))
            {
                ColorMatrix colorMatrix = new(new float[][] { new float[] { .3f, .3f, .3f, 0, 0 }, new float[] { .59f, .59f, .59f, 0, 0 }, new float[] { .11f, .11f, .11f, 0, 0 }, new float[] { 0, 0, 0, 1, 0 }, new float[] { 0, 0, 0, 0, 1 } });

                using ImageAttributes attributes = new();

                attributes.SetColorMatrix(colorMatrix);

                g.DrawImage(original, new System.Drawing.Rectangle(0, 0, original.Width, original.Height), 0, 0, original.Width, original.Height, GraphicsUnit.Pixel, attributes);
            }

            return newBitmap;
        }
    }
}
