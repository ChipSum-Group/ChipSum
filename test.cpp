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

#include "ChipSum.hpp"



int main(int argc, char *argv[]) {

    ChipSum::Common::Init(argc, argv);
    {


    }
    ChipSum::Common::Finalize();
}
