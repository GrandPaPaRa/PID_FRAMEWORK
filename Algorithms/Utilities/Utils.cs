﻿using Emgu.CV;
using Emgu.CV.Structure;
using System;
using System.Collections.Generic;
using System.Linq;
using static System.Math;

namespace Algorithms.Utilities
{
    
    public class Utils
    {
        #region Constants
        public static readonly double B = 0.8;
        #endregion

        #region Change pixel color
        public static void SetPixelColor<TColor>(Image<TColor, byte> inputImage, int row, int column, TColor pixel)
            where TColor : struct, IColor
        {
            if (row >= 0 && row < inputImage.Height && column >= 0 && column < inputImage.Width)
            {
                inputImage[row, column] = pixel;
            }
        }
        #endregion

        #region Merge two images
        public static Image<Bgr, byte> Merge(IImage leftImage, IImage rightImage, int borderWidth = 0)
        {
            Image<Bgr, byte> img1 = ConvertToBgr(leftImage);
            Image<Bgr, byte> img2 = ConvertToBgr(rightImage);

            int maxHeight = Max(img1.Height, img2.Height);
            int maxWidth = Max(img1.Width, img2.Width);

            Image<Bgr, byte> result = new Image<Bgr, byte>(2 * maxWidth + borderWidth, maxHeight);

            int verticalOffset = 0, horizontalOffset = 0;

            if (img1.Height != maxHeight || img1.Width != maxWidth)
            {
                verticalOffset = (maxHeight - img1.Height) / 2;
                horizontalOffset = (maxWidth - img1.Width) / 2;
            }

            for (int y = verticalOffset; y < img1.Height + verticalOffset; ++y)
            {
                for (int x = horizontalOffset; x < img1.Width + horizontalOffset; ++x)
                {
                    result[y, x] = img1[y - verticalOffset, x - horizontalOffset];
                }
            }

            verticalOffset = horizontalOffset = 0;

            if (img2.Height != maxHeight || img2.Width != maxWidth)
            {
                verticalOffset = (maxHeight - img2.Height) / 2;
                horizontalOffset = (maxWidth - img2.Width) / 2;
            }

            for (int y = verticalOffset; y < img2.Height + verticalOffset; ++y)
            {
                for (int x = horizontalOffset + maxWidth + borderWidth; x < img2.Width + horizontalOffset + maxWidth + borderWidth; ++x)
                {
                    result[y, x] = img2[y - verticalOffset, x - horizontalOffset - maxWidth - borderWidth];
                }
            }

            return result;
        }

        private static Image<Bgr, byte> ConvertToBgr(IImage image)
        {
            return image is Image<Gray, byte> grayImg
                ? grayImg.Convert<Bgr, byte>()
                : image as Image<Bgr, byte>;
        }
        #endregion

        #region Compute histogram
        public static int[] ComputeHistogram(Image<Gray, byte> inputImage)
        {
            int[] histogram = new int[256];

            for (int y = 0; y < inputImage.Height; y++)
            {
                for (int x = 0; x < inputImage.Width; x++)
                {
                    ++histogram[inputImage.Data[y, x, 0]];
                }
            }

            return histogram;
        }
        #endregion

        #region Swap
        public static void Swap<T>(ref T lhs, ref T rhs)
        {
            (rhs, lhs) = (lhs, rhs);
        }
        #endregion

        #region Calculate Integral Image 
        public static Image<Gray, double> CalculateIntegralImage(Image<Gray, byte> initialImage) {
            Image<Gray, double> integralImage = new Image<Gray, double>(initialImage.Size);
            for (int y = 0; y < initialImage.Height; y++) {
                for (int x = 0; x < initialImage.Width; x++) {
                    if (x == 0 && y == 0)
                    {
                        integralImage.Data[y, x, 0] = initialImage.Data[y, x, 0];
                    }
                    else if (x == 0)
                    {
                        integralImage.Data[y, x, 0] = integralImage.Data[y - 1, 0, 0] + initialImage.Data[y, 0, 0];
                    }
                    else if (y == 0)
                    {
                        integralImage.Data[y, x, 0] = integralImage.Data[0, x - 1, 0] + initialImage.Data[0, x, 0];
                    }
                    else {
                        integralImage.Data[y, x, 0] = integralImage.Data[y, x - 1, 0] + integralImage.Data[y - 1, x, 0]
                            - integralImage.Data[y - 1, x - 1, 0] + initialImage.Data[y, x, 0];
                    }
                    
                }
            }
            return integralImage;
        }
        #endregion

        #region Median Filter for one pixel
        public static void ApplyMedianFilterAtPixel<TColor>(Image<TColor, byte> inputImage, Image<TColor, byte> processedImage, int x, int y, int windowSize) where TColor : struct, IColor
        {
            // Ensure the window size is odd and greater than 1
            if (windowSize < 3)
            {
                throw new ArgumentException("Window size must be greater than 3.");
            }

            int halfWindowSize = windowSize / 2;

            // Lists to hold values for each channel
            List<byte> redValues = new List<byte>();
            List<byte> greenValues = new List<byte>();
            List<byte> blueValues = new List<byte>();

            // Collect pixel values in the neighborhood
            for (int dy = -halfWindowSize; dy <= halfWindowSize; dy++)
            {
                for (int dx = -halfWindowSize; dx <= halfWindowSize; dx++)
                {
                    int neighborX = x + dx;
                    int neighborY = y + dy;

                    // Check bounds to avoid accessing outside the image
                    if (neighborX >= 0 && neighborX < inputImage.Width && neighborY >= 0 && neighborY < inputImage.Height)
                    {
                        if (typeof(TColor) == typeof(Gray))
                        {
                            // For grayscale images, add intensity directly
                            redValues.Add(inputImage.Data[neighborY, neighborX, 0]);  // Only one channel for grayscale
                        }
                        else if (typeof(TColor) == typeof(Rgb) || typeof(TColor) == typeof(Bgr) )
                        {
                            // For RGB images, add each channel to its respective list
                            redValues.Add(inputImage.Data[neighborY, neighborX, 0]);
                            greenValues.Add(inputImage.Data[neighborY, neighborX, 1]);
                            blueValues.Add(inputImage.Data[neighborY, neighborX, 2]);
                        }
                    }
                }
            }
            byte medianRed = 0, medianGreen = 0, medianBlue = 0;
            // Calculate the median for each channel
            if (typeof(TColor) == typeof(Gray))
            {
                medianRed = CalculateMedian(redValues);
            }
            else {
                medianRed = CalculateMedian(redValues);
                medianGreen = CalculateMedian(greenValues);
                medianBlue = CalculateMedian(blueValues);
            }
            

            // Set the median value back to the pixel
            if (typeof(TColor) == typeof(Gray))
            {
                processedImage.Data[y, x, 0] = medianRed;  // Using redValues since grayscale uses a single list
            }
            else if (typeof(TColor) == typeof(Rgb) || typeof(TColor) == typeof(Bgr))
            {
                processedImage.Data[y, x, 0] = medianRed;
                processedImage.Data[y, x, 1] = medianGreen;
                processedImage.Data[y, x, 2] = medianBlue;

            }
        }

        private static byte CalculateMedian(List<byte> values)
        {
            values.Sort();
            int count = values.Count;
            if (count % 2 == 1)
            {
                // Odd number of elements, take the middle one
                return values[count / 2];
            }
            else
            {
                // Even number of elements, average the two middle ones
                return (byte)((values[(count / 2) - 1] + values[count / 2]) / 2);
            }
        }

        #endregion

        #region Dilation
        public static Image<Bgr, byte> Dilation(Image<Bgr, byte> initialImage, int maskSize) {
            Image<Bgr, byte> result = new Image<Bgr, byte>(initialImage.Size);
            for (int y = 0; y < initialImage.Height; y++) {
                for (int x = 0; x < initialImage.Width; x++) {
                    double maxColorValue = -1;
                    Bgr maxColor = new Bgr(0,0,0);
                    for (int j = -maskSize / 2; j <= maskSize / 2; j++) {
                        for (int i = -maskSize / 2; i <= maskSize / 2; i++)
                        {
                            if (y + i >= 0 && y + i < initialImage.Height && x + j >= 0 && x + j < initialImage.Width) {
                                Bgr pixel = initialImage[y + i, x + j];
                                double currentColorValue = Math.Sqrt(pixel.Red * pixel.Red + pixel.Green * pixel.Green + pixel.Blue * pixel.Blue);
                                if (currentColorValue > maxColorValue) {
                                    maxColor = pixel;
                                    maxColorValue = currentColorValue;
                                }
                            }
                            
                        }
                    }
                    result[y, x] = maxColor;
                }
            }

            return result;
        }
        #endregion
        #region Erosion
        public static Image<Bgr, byte> Erosion(Image<Bgr, byte> initialImage, int maskSize)
        {
            Image<Bgr, byte> result = new Image<Bgr, byte>(initialImage.Size);
            for (int y = 0; y < initialImage.Height; y++)
            {
                for (int x = 0; x < initialImage.Width; x++)
                {
                    double minColorValue = 999999;
                    Bgr minColor = new Bgr(0, 0, 0);
                    for (int j = -maskSize / 2; j <= maskSize / 2; j++)
                    {
                        for (int i = -maskSize / 2; i <= maskSize / 2; i++)
                        {
                            if (y + i >= 0 && y + i < initialImage.Height && x + j >= 0 && x + j < initialImage.Width)
                            {
                                Bgr pixel = initialImage[y + i, x + j];
                                double currentColorValue = Math.Sqrt(pixel.Red * pixel.Red + pixel.Green * pixel.Green + pixel.Blue * pixel.Blue);
                                if (currentColorValue < minColorValue)
                                {
                                    minColor = pixel;
                                    minColorValue = currentColorValue;
                                }
                            }

                        }
                    }
                    result[y, x] = minColor;
                }
            }

            return result;
        }
        #endregion
    }
}