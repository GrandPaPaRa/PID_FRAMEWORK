using Emgu.CV;
using Emgu.CV.Structure;
using System;
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

    }
}