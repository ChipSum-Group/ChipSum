/* * * * * * * * * * * * * * * * * * * * *
 *   File:     cg.cpp
 *   Author:   Yaojie Yu
 *   group:    CDCS-HPC
 *   Time:     2022-3-16
 * * * * * * * * * * * * * * * * * * * * * */

#include "../ChipSum.hpp"
#include "../chipsum/chipsum_macro.h"

Vector gmres(CSR &A, Vector &b, Vector &x, double tol, int max_it)
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

    int m = max_it;
    int n = b.GetSize();

    // get r = b - A*x
    Vector r(n);
    A.SPMV(x, r);
    b.AXPBY(r, 1.0, -1.0); // r = b - A*x

    double bnrm2 = b.Norm2();
    if (bnrm2 == 0.0)
        bnrm2 = 1.0;

    double error = r.Norm2() / bnrm2;

    if (error < tol)
        return x;

    Vector sn(m);
    Vector cs(m);
    Vector beta(m + 1);

    beta(0) = r.Norm2();

    Matrix H(m + 1, m);
    Matrix Q(n, m + 1);

    r *= 1.0 / beta(0);

    Q.SetCol(0, r);

    int j = 0;

    Vector tmp1(n);

    Vector tmp2(n);

    while (j < m)
    {
        // Arnoldi process
        Q.GetColCopy(j, tmp1);
        A.SPMV(tmp1, tmp2);
        Q.SetCol(j+1, tmp2);

        for (int i = 0; i <= j; i++) {
            Q.GetColCopy(j,   tmp1);
            Q.GetColCopy(j+1, tmp2);

            H(i, j) = tmp1.Dot(tmp2);

            Q.GetColCopy(j,   tmp1);
            Q.GetColCopy(j+1, tmp2);

            tmp1 *= -1.0*H(i,j);
            tmp2 += tmp1;

            Q.SetCol(j+1, tmp2);
        }

        Q.GetColCopy(j+1, tmp1);
        H(j+1, j) = tmp1.Norm2();
        
        tmp1 *= 1.0/H(j+1,j);
        Q.SetCol(j+1, tmp1);

        // Applying Givens Rotation to H col
        for (int i = 0; i <= j - 1; i++) {
            double temp = cs(i) * H(i, j) + sn(i) * H(i + 1, j);
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

        if (error <= tol) break;

        j++;
    }


    return x;
}

int main(int argc, char *argv[])
{

    /* .mtx格式数据，如$HOME/ChipSum/data/A.mtx */
    char *filename_A = argv[1];

    /* .mtx格式数据，如$HOME/ChipSum/data/b.csv */
    char *filename_b = argv[2];

    ChipSum::Common::Init(argc, argv);
    {

        CSInt nv = 0, ne = 0;
        CSInt *xadj, *adj;
        double *ew;

        KokkosKernels::Impl::read_matrix<CSInt, CSInt, double>(&nv, &ne, &xadj, &adj, &ew, filename_A);

        CSR A(nv, nv, ne, xadj, adj, ew);

        vector<double> b_data;
        double temp;

        ifstream IN(filename_b);

        for (int i = 0; i < nv; ++i)
        {
            temp = 1.;
            b_data.push_back(temp);
        }

        IN.close();

        Vector b(nv, b_data.data());

        Vector x0(nv);

        x0 *= 0;

        double tol = 1e-12;
        int max_it = 500;

        auto sol_gmres = gmres(A, b, x0, tol, max_it);

        sol_gmres.Print();

        delete xadj;
        delete adj;
        delete ew;
    }
    ChipSum::Common::Finalize();
}