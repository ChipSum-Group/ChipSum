/*
 * @Description: 
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-17 09:59:56
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-17 10:14:22
 */
#ifndef __CHIPSUM_BMP_WRITER_H__
#define __CHIPSUM_BMP_WRITER_H__

namespace ChipSum {
namespace Common {
void WriteBMP(int w,int h,char *img, const char *filename);
} // namespace Common
} // namespace ChipSum

#endif