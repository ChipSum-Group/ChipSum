///
/// \file     spilu.cpp
/// \author   Zhou Xingbin
/// \group    CDCS-HPC
/// \date     2022-3-17
/// \brief    %stuff%
///

#include <iostream>
using namespace std;

#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosKernels_default_types.hpp>

#include "../ChipSum.hpp"
#include "../chipsum/chipsum_macro.h"


int main(int argc, char *argv[]) {

    ChipSum::Common::Init(argc, argv);
    {

        CSInt nrows = 5;
        CSInt nnz   = 10;
        CSInt *row_map = (CSInt *)malloc((nrows + 1) * sizeof(CSInt));
        CSInt *entries = (CSInt *)malloc(nnz * sizeof(CSInt));
        CSFloat *values = (CSFloat *)malloc(nnz * sizeof(CSFloat));
        CSVector b(nrows,1.0);
        CSVector x(nrows);

        // Fill rowmap, entries, of Crs graph for simple example matrix, values set to ones
        //  [ [1 0 1 0 0]
        //    [0 1 0 0 1]
        //    [0 0 1 1 1]
        //    [0 0 0 1 1]
        //    [0 0 0 0 1] ];

        row_map[0] = 0;
        row_map[1] = 2;
        row_map[2] = 4;
        row_map[3] = 7;
        row_map[4] = 9;
        row_map[5] = 10;

        entries[0]  = 0;
        entries[1]  = 2;
        entries[2]  = 1;
        entries[3]  = 4;
        entries[4]  = 2;
        entries[5]  = 3;
        entries[6]  = 4;
        entries[7]  = 3;
        entries[8]  = 4;
        entries[9]  = 4;

        values[0]  = 1;
        values[1]  = 1;
        values[2]  = 1;
        values[3]  = 1;
        values[4]  = 1;
        values[5]  = 1;
        values[6]  = 1;
        values[7]  = 1;
        values[8]  = 1;
        values[9]  = 1;
        
        CSR A(nrows,nrows,nnz,row_map,entries,values);
        

        A.SPTRSV(x,b,false);
        std::cout<<"Upper triangular solve:"<<std::endl;
        x.Print();

        
        // Fill rowmap, entries, of Crs graph for simple example matrix, values set to ones

       //  [ [1 0 0 0 0]
       //    [0 1 0 0 0]
       //    [1 0 1 0 0]
       //    [0 0 1 1 0]
       //    [1 1 1 1 1] ];

        row_map[0] = 0;
        row_map[1] = 1;
        row_map[2] = 2;
        row_map[3] = 4;
        row_map[4] = 6;
        row_map[5] = 10;
        
        entries[0] = 0;
        entries[1] = 1;
        entries[2] = 0;
        entries[3] = 2;
        entries[4] = 2;
        entries[5] = 3;
        entries[6] = 1;
        entries[7] = 2;
        entries[8] = 3;
        entries[9] = 4;
        CSR A2(nrows,nrows,nnz,row_map,entries,values);
        A2.SPTRSV(x,b,true);
        std::cout<<"Lower triangular solve:"<<std::endl;
        x.Print();
        
        free(row_map);
        free(entries);
        free(values);
    }
    ChipSum::Common::Finalize();
}
