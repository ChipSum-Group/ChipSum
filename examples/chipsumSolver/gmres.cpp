/* * * * * * * * * * * * * * * * * * * * *
 *   File:     cg.cpp
 *   Author:   Yaojie Yu
 *   group:    CDCS-HPC
 *   Time:     2022-06-07
 * * * * * * * * * * * * * * * * * * * * * */

#include <gmres.hpp>

int main(int argc, char *argv[])
{

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
        CSFloat *values = (CSFloat *)malloc(13 * sizeof(CSFloat));

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

        CSR A(nrows, ncols, annz, row_map, col_map, values);

        // A.PrintPattern();

        CSFloat *vec1 = new CSFloat[5];
        CSFloat *vec2 = new CSFloat[5];

        for (int i = 0; i < 5; ++i)
        {
            vec1[i] = i * i + 1.0; // double(i);
            vec2[i] = 0.0;
        }

        CSVector b(5, vec1);
        CSVector x0(5, vec2);

        cout << "A = ";
        A.Print();

        cout << "b = ";
        b.Print();

        cout << "x0 = ";
        x0.Print();

        CSFloat tol = 1e-12;
        CSInt max_it = 20;

        auto sol_gmres = ChipSum::Solver::gmres(A, b, x0, tol, max_it);

        cout << "sol_gmres = ";
        sol_gmres.Print();
        cout << " " << endl;

        CSVector res(n);
        cout << "A * sol_gmres= ";
        A.SPMV(sol_gmres, res);
        res.Print();
        cout << " " << endl;

        std::free(row_map);
        std::free(col_map);
        std::free(values);

        delete[] vec1;
        delete[] vec2;
    }

    ChipSum::Common::Finalize();
}