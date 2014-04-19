


#include <iostream>

#include "TiffUtilities.h"


/**
 * @brief This program takes a single argument, the path to the 16 bit image
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cout << "This program needs a single argument which is the path to the tiff iamge" << std::endl;
    return EXIT_FAILURE;
  }

  // Allocate a TiffImage Structure
  TiffImage* tiffImage = static_cast<TiffImage*>(malloc(sizeof(TiffImage)));

  // Create an instance of the TiffUtilities Class
  TiffUtilities tifUtil;

  int err = tifUtil.readInputImage(tiffImage, argv[1]);
  if (err < 0)
  {
    return EXIT_FAILURE;
  }

  // At this point the image data is stored in the tiffImage->imageData pointer
  // which is a void pointer. If you are going to work with it you need to know
  // how to cast it. This next if-then-else does the decision making for you. THere
  // are other ways of doing this decision making process. C++ with Templates might
  // be a slightly more efficient way of doing this.


  if (tiffImage->bitsPerSample == 16)
  {
    std::cout << "16 Bit Tiff Image" << std::endl;
    unsigned short max = 0;
    unsigned short min = 0xFFFF;
    unsigned short* pixels = static_cast<unsigned short*>(tiffImage->imageData);
    size_t totalPixels = tiffImage->width * tiffImage->height;
    for(size_t p = 0; p < totalPixels; ++p)
    {
      if (pixels[p] > max) { max = pixels[p]; };
      if (pixels[p] < min) { min = pixels[p]; };
    }
    std::cout << "Max 16 Bit Value: " << max << std::endl;
    std::cout << "Min 16 Bit Value: " << min << std::endl;
  }
  else if(tiffImage->bitsPerSample == 8)
  {
    std::cout << "8 Bit Tiff Image" << std::endl;
  }
  else
  {
    std::cout << "Unknown Tiff bitsPerSample" << std::endl;
  }




  // When you are completely done with the image data you need to free it:
  free(tiffImage->imageData);

  return EXIT_SUCCESS;
}


