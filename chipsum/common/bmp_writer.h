///
/// \file     bmp_writer.h
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-11-01
/// \brief    BMP图片功能类，目前我们用这两个接口实现对线性系统
///           pattern的可视化输出.
///
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
