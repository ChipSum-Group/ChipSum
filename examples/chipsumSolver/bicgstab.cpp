/* * * * * * * * * * * * * * * * * * * * *
 *   File:     bicg.cpp
 *   Author:   Yaojie Yu
 *   group:    CDCS-HPC
 *   Time:     2021-12-14
 * * * * * * * * * * * * * * * * * * * * * */

#include "../ChipSum.hpp"
#include "../chipsum/chipsum_macro.h"

CSVector bicgstab(CSR &A, CSVector &b, CSVector &x, CSFloat tol, int max_it)
{
    // BiConjugate Gradient Stabilized Method without preconditioning.
    //
    // input   A        REAL matrix
    //         x        REAL initial guess vector
    //         b        REAL right hand side vector
    //         tol      REAL error tolerance
    //         max_it   INTEGER maximum number of iterations
    //
    // output  x        REAL solution vector

    CSFloat bnrm2 = b.Norm2();
    if (bnrm2 == 0.0)
        bnrm2 = 1.0;

    CSVector r(x.GetSize());
    A.SPMV(x, r);
    b.AXPBY(r, 1.0, -1.0); // r = b - A*x

    CSFloat error = r.Norm2() / bnrm2;

    if (error < tol)
        return x;

    CSVector r_tld(x.GetSize());
    r_tld.DeepCopy(r);

    CSVector p(x.GetSize()), p_hat(x.GetSize());
    CSVector v(x.GetSize()), t(x.GetSize());
    CSVector s(x.GetSize()), s_hat(x.GetSize());

    CSFloat omega = 1.0;
    CSFloat alpha, beta, rho, rho_1;

    for (int i = 0; i < max_it; i++)
    {
        rho = r.Dot(r_tld);

        if (rho == 0.0)
            break;

        if (i > 0)
        {
            beta = (rho / rho_1) * (alpha / omega);
            v.AXPBY(p, -omega, 1.0); // p = p - omega* v
            r.AXPBY(p, 1.0, beta);   //   p = r + beta * p
        }
        else
        {
            p.DeepCopy(r); // p = r
        }

        p_hat.DeepCopy(p); // p_hat = M^-1 * p, M is a preconditioner matrix
        A.SPMV(p_hat, v);  // v = A* p_hat
        alpha = rho / (r_tld.Dot(v));

        s.DeepCopy(v);           // s = v
        r.AXPBY(s, 1.0, -alpha); // s = r - alpha * v

        if (s.Norm2() < tol)
        {
            p_hat.AXPBY(x, alpha, 1.0); // x = x + alpha * p_hat
            break;
        }

        s_hat.DeepCopy(s); // s_hat = M^-1 * s, M is a preconditioner matrix

        A.SPMV(s_hat, t); // t = A*s_hat

        CSFloat tmp = t.Dot(t);
        omega = (t.Dot(s)) / tmp;
        // omega = (t.Dot(s)) / (t.Dot(t));
        p_hat.AXPBY(x, alpha, 1.0); // x = x +  alpha * p_hat
        s_hat.AXPBY(x, omega, 1.0); // x = x +  omega * s_hat

        r.DeepCopy(s);
        t.AXPBY(r, -omega, 1.0); // r = s - omega * t

        error = r.Norm2() / bnrm2;

        if (error <= tol)
            break;
        if (omega == 0.0)
            break;

        rho_1 = rho;
        CSVector r_temp(b.GetSize());
        A.SPMV(x, r_temp);          /* r_temp = A*x */
        b.AXPBY(r_temp, 1.0, -1.0); /* r_temp = b-r_temp */
    }

    return x;
}

// int main(int argc, char *argv[])
// {

//     char *filename_A = argv[1];
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
//         int max_it = 500;

//         auto sol_bicgstab = bicgstab(A, b, x0, tol, max_it);

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

        auto sol = bicgstab(A, b, x0, tol, max_it);

        cout << "sol = ";
        sol.Print();
        cout << " " << endl;

        CSVector res(n);
        cout << "A * sol= ";
        A.SPMV(sol, res);
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