#ifndef __CHIPSUM_MACRO_H__
#define __CHIPSUM_MACRO_H__



#define CHIPSUM_FUNCTION_INLINE inline /**<---根据需求修改*/
#define CHIPSUM_DECLARED_FUNCTION inline /**<---根据需求修改*/


#define CHIPSUM_UNUSED(x) (void)x;


#define CSERR_t int

#define CSInt_t int

//#ifdef USE_CUDA
//#define CHIPSUM_FUNCTION_INLINE __global__ __device__
//#endif

#endif

