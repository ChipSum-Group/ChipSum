/*
 * @Description: 
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-17 10:00:16
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-17 10:14:33
 */

#include <cstdio>
#include <cstdlib>

#include "bmp_writer.h"


void ChipSum::Common::WriteBMP(int w,int h,char *img, const char *filename) {
  int l = (w * 3 + 3) / 4 * 4;
  int bmi[] = {1 * h + 54, 0,     54, 40, w,   h, 1 | 3 * 8 << 16,
               0,          l * h, 0,  0,  100, 0};
  FILE *fp = std::fopen(filename, "wb");
  std::fprintf(fp, "BM");
  std::fwrite(&bmi, 52, 1, fp);
  std::fwrite(img, 1, l * h, fp);
  std::fclose(fp);
}