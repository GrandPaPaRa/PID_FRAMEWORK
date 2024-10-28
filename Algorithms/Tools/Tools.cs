using Emgu.CV;
using Emgu.CV.Structure;
using System.Windows.Forms;
using System.Windows;
using System.Drawing;
using System;
namespace Algorithms.Tools
{
    public class Tools
    {
        #region Copy
        public static Image<Gray, byte> Copy(Image<Gray, byte> inputImage)
        {
            Image<Gray, byte> result = inputImage.Clone();
            return result;
        }

        public static Image<Bgr, byte> Copy(Image<Bgr, byte> inputImage)
        {
            Image<Bgr, byte> result = inputImage.Clone();
            return result;
        }
        #endregion

        #region Invert
        public static Image<Gray, byte> Invert(Image<Gray, byte> inputImage)
        {
            Image<Gray, byte> result = new Image<Gray, byte>(inputImage.Size);
            for (int y = 0; y < inputImage.Height; y++)
            {
                for (int x = 0; x < inputImage.Width; x++)
                {
                    result.Data[y, x, 0] = (byte)(255 - inputImage.Data[y, x, 0]);
                }
            }
            return result;
        }

        public static Image<Bgr, byte> Invert(Image<Bgr, byte> inputImage)
        {
            Image<Bgr, byte> result = new Image<Bgr, byte>(inputImage.Size);
            for (int y = 0; y < inputImage.Height; y++)
            {
                for (int x = 0; x < inputImage.Width; x++)
                {
                    result.Data[y, x, 0] = (byte)(255 - inputImage.Data[y, x, 0]);
                    result.Data[y, x, 1] = (byte)(255 - inputImage.Data[y, x, 1]);
                    result.Data[y, x, 2] = (byte)(255 - inputImage.Data[y, x, 2]);
                }
            }
            return result;
        }
        #endregion

        #region Convert color image to grayscale image
        public static Image<Gray, byte> Convert(Image<Bgr, byte> inputImage)
        {
            Image<Gray, byte> result = inputImage.Convert<Gray, byte>();
            return result;
        }
        #endregion

        #region Binary
        public static Image<Gray, byte> Binary(Image<Gray, byte> inputImage, double threshold) { 
            Image<Gray, byte> result = new Image<Gray, byte>(inputImage.Size);
            for (int y = 0; y<inputImage.Height; ++y) { 
                for (int x = 0; x<inputImage.Width; ++x) { 
                    if (inputImage.Data[y, x, 0] > threshold) result.Data[y, x, 0] = 255; 
                }
            }
            return result;
        }
        #endregion

        #region Crop
        public static Image<TColor, TDepth> Crop<TColor, TDepth>(
        Image<TColor, TDepth> inputImage,
        Point firstMousePos,
        Point lastMousePos)
        where TColor : struct, IColor
        where TDepth : struct
        {
            // Determine leftTop and rightBottom points based on mouse positions
            Point leftTop = new Point(
                Math.Min(firstMousePos.X, lastMousePos.X),
                Math.Min(firstMousePos.Y, lastMousePos.Y)
            );
            Point rightBottom = new Point(
                Math.Max(firstMousePos.X, lastMousePos.X),
                Math.Max(firstMousePos.Y, lastMousePos.Y)
            );

            // Calculate crop width and height
            int cropWidth = rightBottom.X - leftTop.X;
            int cropHeight = rightBottom.Y - leftTop.Y;

            // Validate crop area dimensions
            if (cropWidth <= 0 || cropHeight <= 0 ||
                leftTop.X < 0 || leftTop.Y < 0 ||
                rightBottom.X > inputImage.Width || rightBottom.Y > inputImage.Height)
            {
                throw new ArgumentException("Invalid crop boundaries.");
            }

            // Set the ROI
            inputImage.ROI = new Rectangle(leftTop.X, leftTop.Y, cropWidth, cropHeight);

            
            
            // Create a new instance of the cropped image
            Image<TColor, TDepth> result = new Image<TColor, TDepth>(cropWidth, cropHeight);
            inputImage.CopyTo(result);  // Copy the cropped region to the result

            // Reset the ROI
            inputImage.ROI = Rectangle.Empty;

            //Display mean and rms
            (double mean, double rms) = CalculateMeanAndRMS(inputImage,
                new Rectangle(leftTop.X, leftTop.Y, cropWidth, cropHeight));

            MessageBox.Show($"Media: {mean}\nAbaterea Medie Pătratică: {rms}", "Rezultate");
            return result;
        }
        #endregion

        #region CalculateMeanAndRMS
        public static (double mean, double rms) CalculateMeanAndRMS<TColor, TDepth>(
       Image<TColor, TDepth> inputImage, Rectangle cropArea)
       where TColor : struct, IColor
       where TDepth : struct
        {
            // Set the ROI (Region of Interest) to the specified crop area
            inputImage.ROI = cropArea;

            double sum = 0;
            double sumSquared = 0;
            int pixelCount = cropArea.Width * cropArea.Height;

            // Iterate through the pixels in the cropped area
            for (int y = 0; y < cropArea.Height; y++)
            {
                for (int x = 0; x < cropArea.Width; x++)
                {
                    // Get pixel value
                    TColor pixelValue = inputImage[y, x];

                    // Calculate intensity based on color type
                    double intensity;
                    if (typeof(TColor) == typeof(Bgr))
                    {
                        var bgr = (Bgr)(object)pixelValue; // Cast to Bgr
                        intensity = 0.114 * bgr.Blue + 0.587 * bgr.Green + 0.299 * bgr.Red; // Calculate intensity
                    }
                    else if (typeof(TColor) == typeof(Rgb))
                    {
                        var rgb = (Rgb)(object)pixelValue; // Cast to Rgb
                        intensity = 0.299 * rgb.Red + 0.587 * rgb.Green + 0.114 * rgb.Blue; // Calculate intensity
                    }
                    else if (typeof(TColor) == typeof(Gray))
                    {
                        var gray = (Gray)(object)pixelValue; // Cast to Gray
                        intensity = gray.Intensity; // Get intensity directly for Gray
                    }
                    else
                    {
                        throw new InvalidOperationException("Unsupported color type");
                    }

                    sum += intensity; // Accumulate the sum
                    sumSquared += intensity * intensity; // Accumulate the squared sum
                }
            }

            // Calculate the mean
            double mean = sum / pixelCount;

            // Calculate the RMS
            double meanSquare = sumSquared / pixelCount;
            double rms = Math.Sqrt(meanSquare - (mean * mean));

            // Reset ROI
            inputImage.ROI = Rectangle.Empty;

            return (mean, rms);
        }
        #endregion
    }
    
}