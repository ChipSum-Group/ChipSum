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

        CSInt nrows = 4;
        CSInt nnz   = 8;
        CSInt *row_map = (CSInt *)malloc((nrows + 1) * sizeof(CSInt));
        CSInt *entries = (CSInt *)malloc(nnz * sizeof(CSInt));
        CSFloat *values = (CSFloat *)malloc(nnz * sizeof(CSFloat));
        CSInt *entries2 = (CSInt *)malloc(nnz * sizeof(CSInt));
        CSFloat *values2 = (CSFloat *)malloc(nnz * sizeof(CSFloat));

        row_map[0] = 0;
        row_map[1] = 2;
        row_map[2] = 4;
        row_map[3] = 6;
        row_map[4] = 8;

        entries[0]  = 0;
        entries[1]  = 2;
        entries[2]  = 0;
        entries[3]  = 2;
        entries[4]  = 0;
        entries[5]  = 2;
        entries[6]  = 0;
        entries[7]  = 2;
        

        values[0]  = 1;
        values[1]  = 1;
        values[2]  = 1;
        values[3]  = 1;
        values[4]  = 1;
        values[5]  = 1;
        values[6]  = 1;
        values[7]  = 1;
        

        entries2[0]  = 1;
        entries2[1]  = 3;
        entries2[2]  = 1;
        entries2[3]  = 3;
        entries2[4]  = 1;
        entries2[5]  = 3;
        entries2[6]  = 1;
        entries2[7]  = 3;
        

        values2[0]  = 2;
        values2[1]  = 2;
        values2[2]  = 2;
        values2[3]  = 2;
        values2[4]  = 2;
        values2[5]  = 2;
        values2[6]  = 2;
        values2[7]  = 2;
        
        
        CSR A(nrows,nrows,nnz,row_map,entries,values);
        CSR B(nrows,nrows,nnz,row_map,entries2,values2);
        CSR C;

        A.Print();
        A.PrintPattern();
        B.Print();
        B.PrintPattern();
        
        // C.SPADD(A,B) = C.SPADD(1.0,A,1.0,B) = C.SPADD(1.0,A,1.0,B,false)
        // C.SPADD(2.0,A,1.0,B) = C.SPADD(2.0,A,1.0,B)
        C.SPADD(2.0,A,1.0,B,true);
        
        C.Print();
        C.PrintPattern();
              
        
        free(row_map);
        free(entries);
        free(entries2);
        free(values);
        free(values2);
    }
    ChipSum::Common::Finalize();
}
