#ifndef __CHIPSUM_MACRO_H__
#define __CHIPSUM_MACRO_H__

#include <cctype>


/// \brief 根据需求修改
#define CHIPSUM_FUNCTION_INLINE inline

/// \brief 根据需求修改
#define CHIPSUM_DECLARED_FUNCTION inline

#if defined(ChipSum_USE_KokkosKernels) || defined(ChipSum_USE_KokkosKernels64)
#define CHIPSUM_SPECIAL_INLINE KOKKOS_INLINE_FUNCTION
// #elif
#else
#define CHIPSUM_SPECIAL_INLINE inline
#endif
#define CHIPSUM_UNUSED(x) (void)x;

#define CSErr_t __int32_t


/// \brief 默认的宏定义类型
#define CSInt __int32_t
#define CSFloat double



#endif

