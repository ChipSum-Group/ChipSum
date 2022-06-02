/* * * * * * * * * * * * * * * * * * * * *
 *   File:     cg.cpp
 *   Author:   Yaojie Yu
 *   group:    CDCS-HPC
 *   Time:     2022-3-16
 * * * * * * * * * * * * * * * * * * * * * */

#include "../ChipSum.hpp"
#include "../chipsum/chipsum_macro.h"

CSVector gmres(CSR &A, CSVector &b, CSVector &x, CSFloat tol, int max_it)
{
    // GMRES Method without preconditioning.
    //
    // input   A        REAL matrix
    //         x        REAL initial guess vector
    //         b        REAL right hand side vector
    //         tol      REAL error tolerance
    //         max_it   INTEGER maximum number of iterations
    //
    // output  x        REAL solution vector

    // get norm of b
    CSFloat bnrm2 = b.Norm2();

    if (bnrm2 == 0.0) bnrm2 = 1.0;

    int n = b.GetSize();
    int m = max_it;

    CSVector r(n);
    A.SPMV(x, r);
    b.AXPBY(r, 1.0, -1.0); // r = b - A*x

    CSFloat error = r.Norm2() / bnrm2;

    if (error < tol) return x;

    CSVector sn(m);
    CSVector cs(m);

    CSVector beta(m + 1);
    beta(0) = r.Norm2();

    CSMatrix H(m + 1, m);
    CSMatrix Q(n, m + 1);

    r *= (1.0 / r.Norm2());
    Q.SetCol(0, r);

    CSVector tmp1(n);
    CSVector tmp2(n);
    CSVector beta_n(n);

    int j = 0;
    while (j < m)
    {
        // Arnoldi process
        Q.GetColCopy(j, tmp1);
        A.SPMV(tmp1, tmp2);
        Q.SetCol(j + 1, tmp2);

        for (int i = 0; i <= j; i++)
        {
            Q.GetColCopy(i, tmp1);
            Q.GetColCopy(j + 1, tmp2);

            H(i, j) = tmp1.Dot(tmp2);

            Q.GetColCopy(i, tmp1);
            tmp1 *= (-1.0 * H(i, j));

            Q.GetColCopy(j + 1, tmp2);
            tmp2 += tmp1;

            Q.SetCol(j + 1, tmp2);
        }

        Q.GetColCopy(j + 1, tmp1);
        H(j + 1, j) = tmp1.Norm2();

        tmp1 *= (1.0 / H(j + 1, j));
        Q.SetCol(j + 1, tmp1);

        // Applying Givens Rotation to H col
        for (int i = 0; i <= j - 1; i++)
        {
            CSFloat temp = cs(i) * H(i, j) + sn(i) * H(i + 1, j);
            H(i + 1, j) = -sn(i) * H(i, j) + cs(i) * H(i + 1, j);
            H(i, j) = temp;
        }

        cs(j) = H(j, j) / sqrt(H(j, j) * H(j, j) + H(j + 1, j) * H(j + 1, j));
        sn(j) = H(j + 1, j) / sqrt(H(j, j) * H(j, j) + H(j + 1, j) * H(j + 1, j));

        H(j, j) = cs(j) * H(j, j) + sn(j) * H(j + 1, j);

        H(j + 1, j) = 0.0;

        // update the residual vector
        beta(j + 1) = -sn(j) * beta(j);
        beta(j) = cs(j) * beta(j);

        error = abs(beta(j + 1)) / b.Norm2();

        if (error <= tol)
            break;

        j++;
    }

    // solve the triangular equation
    CSMatrix H_n(j + 1, j + 1);
    H.GetPartSlice(0, 0, j + 1, j + 1, H_n);

    CSMatrix Y(j + 1, 1);

    for (int i = 0; i <= j; i++)
    {
        Y(i, 0) = beta(i);
    }

    Y.TRMM(H_n, 1.0, "L", "U");
    CSVector y(j + 1);
    Y.GetColCopy(0, y);

    CSMatrix Q_n(n, j + 1);
    Q.GetPartSlice(0, 0, n, j + 1, Q_n);

    // get solution x
    CSVector x_n(n);
    Q_n.GEMV(y, x_n);
    x += x_n;

    return x;
}

// int main(int argc, char *argv[])
// {

//     /* .mtx格式数据，如$HOME/ChipSum/data/A.mtx */
//     char *filename_A = argv[1];

//     /* .mtx格式数据，如$HOME/ChipSum/data/b.csv */
//     char *filename_b = argv[2];

//     ChipSum::Common::Init(argc, argv);
//     {

//         CSInt nv = 0, ne = 0;
//         CSInt *xadj, *adj;
//         CSFloat *ew;

//         KokkosKernels::Impl::read_matrix<CSInt, CSInt, CSFloat>(&nv, &ne, &xadj, &adj, &ew, filename_A);

//         CSR A(nv, nv, ne, xadj, adj, ew);

//         vector<CSFloat> b_data;
//         CSFloat temp;

//         ifstream IN(filename_b);

//         for (int i = 0; i < nv; ++i)
//         {
//             temp = 1.;
//             b_data.push_back(temp);
//         }

//         IN.close();

//         CSVector b(nv, b_data.data());

//         CSVector x0(nv);

//         x0 *= 0;

//         CSFloat tol = 1e-12;
//         int max_it = 80;

//         auto sol_gmres = gmres(A, b, x0, tol, max_it);

//         CSVector res(nv);
//         A.SPMV(sol_gmres, res);

//         //res.Print();
//         std::cout << res.Norm2() << std::endl;
//         // sol_gmres.Print();

//         delete xadj;
//         delete adj;
//         delete ew;
//     }
//     ChipSum::Common::Finalize();
// }

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

        auto sol_gmres = gmres(A, b, x0, tol, max_it);

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