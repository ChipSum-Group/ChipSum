/*
 * @Description: 
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-17 09:59:56
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-18 11:38:48
 */
#ifndef __CHIPSUM_BMP_WRITER_H__
#define __CHIPSUM_BMP_WRITER_H__

namespace ChipSum {
namespace Common {

///
/// \brief 输出一个bmp格式的图片，w和h代表图片像素，img储存rgb值，图片名称为filename
/// \param w 图片宽
/// \param h 图片高
/// \param img 图片rgb信息
/// \param filename 图片名称
///
void write_bmp(int w,int h,unsigned char *img, const char *filename);


///
/// \brief 翻转bmp图片
/// \param w 图片宽
/// \param h 图片高
/// \param img 图片rgb信息
///
void flip_bmp(int w,int h,unsigned char *img);
} // namespace Common


} // namespace ChipSum


#endif
