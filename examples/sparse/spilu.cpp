///
/// \file     spilu.cpp
/// \author   Zhou Xingbin
/// \group    CDCS-HPC
/// \date     2022-2-15
/// \brief    %stuff%
///

#include <iostream>
using namespace std;

#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosKernels_default_types.hpp>

#include <vector>

#include "../ChipSum.hpp"
#include "../chipsum/chipsum_macro.h"


int main(int argc, char *argv[]) {

    ChipSum::Common::Init(argc, argv);
    {

        CSInt nrows = 9;
        CSInt nnz   = 21;
        CSInt *row_map = (CSInt *)malloc((nrows + 1) * sizeof(CSInt));
        CSInt *entries = (CSInt *)malloc(nnz * sizeof(CSInt));
        CSFloat *values = (CSFloat *)malloc(nnz * sizeof(CSFloat));

        row_map[0] = 0;
        row_map[1] = 3;
        row_map[2] = 5;
        row_map[3] = 6;
        row_map[4] = 9;
        row_map[5] = 11;
        row_map[6] = 13;
        row_map[7] = 15;
        row_map[8] = 18;
        row_map[9] = nnz;

        entries[0]  = 0;
        entries[1]  = 2;
        entries[2]  = 5;
        entries[3]  = 1;
        entries[4]  = 6;
        entries[5]  = 2;
        entries[6]  = 0;
        entries[7]  = 3;
        entries[8]  = 4;
        entries[9]  = 0;
        entries[10] = 4;
        entries[11] = 1;
        entries[12] = 5;
        entries[13] = 2;
        entries[14] = 6;
        entries[15] = 3;
        entries[16] = 4;
        entries[17] = 7;
        entries[18] = 3;
        entries[19] = 4;
        entries[20] = 8;

        values[0]  = 10;
        values[1]  = 0.3;
        values[2]  = 0.6;
        values[3]  = 11;
        values[4]  = 0.7;
        values[5]  = 12;
        values[6]  = 5;
        values[7]  = 13;
        values[8]  = 1;
        values[9]  = 4;
        values[10] = 14;
        values[11] = 3;
        values[12] = 15;
        values[13] = 7;
        values[14] = 16;
        values[15] = 6;
        values[16] = 5;
        values[17] = 17;
        values[18] = 2;
        values[19] = 2.5;
        values[20] = 18;
        
        CSR A(nrows,nrows,nnz,row_map,entries,values);
        CSR L;
        CSR U;
        A.SPILU(L,U,2);
        std::cout<<"L_CSR:";
        L.Print();
        std::cout<<"U_CSR:";
        U.Print();

        CSVector bb(nrows);
        CSVector bb_tmp(nrows);
        std::vector<CSFloat> one(nrows,1.0);
        CSVector e_one(nrows,one.data());

        A.SPMV(e_one,bb);
        // CSFloat bb_nrm = bb.Norm2();
        // bb.Print();

        U.SPMV(e_one,bb_tmp);
        // CSFloat bb_tmp_nrm = bb_tmp.Norm2();
        // std::cout<<"L_nrm: "<<bb_tmp_nrm<<std::endl; 
        // bb_tmp.Print();

        L.SPMV(bb_tmp,bb,1.0,-1.0);
        // bb.Print();
        CSFloat diff_nrm = bb.Norm2();
        // std::cout<<"bb_nrm: "<< bb_nrm << std::endl;
        std::cout<<"diff_nrm: "<< diff_nrm << std::endl;        
        
        free(row_map);
        free(entries);
        free(values);
    }
    ChipSum::Common::Finalize();
}
