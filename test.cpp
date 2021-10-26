/*
 * @Author       : your name
 * @Date         : 2021-08-10 15:35:49
 * @LastEditTime: 2021-10-11 09:06:48
 * @LastEditors: Li Kunyun
 * @Description  : In User Settings Edit
 * @FilePath     : \\lky\\ChipSum\\test.cpp
 */

#include <iostream>
using namespace std;

// #include <KokkosKernels_IOUtils.hpp>
// #include <KokkosSparse_CrsMatrix.hpp>

#include <type_traits>
#include <vector>

#include "ChipSumConfig.h"
#include "chipsum/backend/backend.hpp"
#include "chipsum/common/enviroment.hpp"
#include "chipsum/numeric/dense_matrix.hpp"

#include "chipsum/numeric/scalar.hpp"
#include "chipsum/numeric/sparse_matrix.hpp"
#include "chipsum/numeric/vector.hpp"



#define N 5

typedef ChipSum::Numeric::Vector<double, std::size_t,
                                 ChipSum::Backend::Kokkos>
    dVector;

typedef ChipSum::Numeric::Scalar<double,std::size_t,ChipSum::Backend::Cuda> tScalar;

int main(int argc, char *argv[]) {

  ChipSum::Common::Init(argc, argv);
  {
      tScalar s(3.2);

      cout<<s()<<endl;
  }
  ChipSum::Common::Finalize();
}
