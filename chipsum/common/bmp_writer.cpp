///
/// \file     bmp_writer.cpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-05
/// \brief    %stuff%
///

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "bmp_writer.h"



void ChipSum::Common::write_bmp(int w, int h,unsigned char *img, const char *filename) {
  int l = (w * 3 + 3) / 4 * 4;
  int bmi[] = {1 * h + 54, 0,     54, 40, w,   h, 1 | 3 * 8 << 16,
               0,          l * h, 0,  0,  100, 0};
  FILE *fp = std::fopen(filename, "wb");
  ::std::fprintf(fp, "BM");
  ::std::fwrite(&bmi, 52, 1, fp);
  ::std::fwrite(img, 1, l * h, fp);

  
  ::std::fclose(fp);
}

void ChipSum::Common::flip_bmp(int w, int h, unsigned char *img) {
  int half = h / 2;
  for (int i = 0; i < half; ++i) {
    for (int j = 0; j < w; ++j) {
      for (int k = 0; k < 3; ++k) {
        ::std::swap(img[i * 3 * w + 3 * j + k],
                  img[(h - i - 1) * 3 * w + 3 * j + k]);
      }
    }
  }
}
