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



#include <type_traits>
#include <vector>




#include "ChipSum.hpp"
#include "chipsum/chipsum_macro.h"




int main(int argc, char *argv[]) {


    ChipSum::Common::Init(argc, argv);
    {


        /*
                *
                *  |  1  0  2  3  0  |
                *  |  0  4  0  5  0  |
                *  |  2  0  6  0  7  |
                *  |  3  5  0  8  0  |
                *  |  0  0  7  0  9  |
                *
                *  nrows = 5
                *  ncols = 5
                *  annz = 13
                *  row_map = {0,3,5,8,11,13}
                *  col_map = {0,2,3,1,3,0,2,4,0,1,3,2,4}
                *  values = {1,2,3,4,5,2,6,7,3,5,8,7,9}
                *
                *  x = {1,1,1,1,1}
                *  y = {6,9,15,16,16}
                *
               */
        CSInt m = 5;
        CSInt n = 5;

        CSInt nrows = m;
        CSInt ncols = n;
        CSInt annz = 13;
        CSInt *row_map = (CSInt *)malloc(6 * sizeof(CSInt));
        CSInt *col_map = (CSInt *)malloc(13 * sizeof(CSInt));
        double *values = (double *)malloc(13 * sizeof(double));
        row_map[0] = 0;
        row_map[1] = 3;
        row_map[2] = 5;
        row_map[3] = 8;
        row_map[4] = 11;
        row_map[5] = 13;
        col_map[0] = col_map[5] = col_map[8] = 0;
        col_map[1] = col_map[6] = col_map[11] = 2;
        col_map[2] = col_map[4] = col_map[10] = 3;
        col_map[3] = col_map[9] = 1;
        col_map[7] = col_map[12] = 4;

        for (int i = 0; i < 5; ++i)
            values[i] = double(i + 1);

        values[5] = 2;
        values[6] = 6;
        values[7] = 7;
        values[8] = 3;
        values[9] = 5;
        values[10] = 8;
        values[11] = 7;
        values[12] = 9;


        //________________________

        CSR s_A(nrows,ncols,annz,row_map,col_map,values);

        s_A.PrintPattern();

        double * vec = new double[5];

        Vector d_vec(5,vec);

        Vector new_vec= s_A*d_vec;

        new_vec.Print();

        std::free(row_map);
        std::free(col_map);
        std::free(values);

        delete [] vec;

    }
    ChipSum::Common::Finalize();
}
