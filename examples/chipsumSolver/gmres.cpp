/* * * * * * * * * * * * * * * * * * * * *
 *   File:     cg.cpp
 *   Author:   Yaojie Yu
 *   group:    CDCS-HPC
 *   Time:     2022-3-16
 * * * * * * * * * * * * * * * * * * * * * */

#include "../ChipSum.hpp"
#include "../chipsum/chipsum_macro.h"

CSVector gmres(CSR &A, CSVector &b, CSVector &x, double tol, int max_it)
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

    int n = b.GetSize();
    int m = max_it;

    // get r = b - A*x
    CSVector r(n);
    A.SPMV(x, r);
    b.AXPBY(r, 1.0, -1.0); // r = b - A*x

    // get norm of b
    double bnrm2 = b.Norm2();

    if (bnrm2 == 0.0)
        bnrm2 = 1.0;

    double error = r.Norm2() / bnrm2;

    if (error < tol)
        return x;

    CSVector sn(m);
    CSVector cs(m);
    CSVector beta(m + 1);

    beta(0) = r.Norm2();

    CSMatrix H(m + 1, m);
    CSMatrix Q(n, m + 1);

    r *= 1.0 / beta(0);

    Q.SetCol(0, r);

    int j = 0;

    CSVector tmp1(n);
    CSVector tmp2(n);
    CSVector beta_n(n);

    for (int j = 0; j < m; j++)
    {
        // Arnoldi process
        Q.GetColCopy(j, tmp1);
        A.SPMV(tmp1, tmp2);
        Q.SetCol(j + 1, tmp2);

        for (int i = 0; i <= j; i++)
        {
            Q.GetColCopy(j, tmp1);
            Q.GetColCopy(j + 1, tmp2);

            H(i, j) = tmp1.Dot(tmp2);

            Q.GetColCopy(j, tmp1);
            Q.GetColCopy(j + 1, tmp2);

            tmp1 *= -1.0 * H(i, j);
            tmp2 += tmp1;

            Q.SetCol(j + 1, tmp2);
        }

        Q.GetColCopy(j + 1, tmp1);
        H(j + 1, j) = tmp1.Norm2();

        tmp1 *= 1.0 / H(j + 1, j);
        Q.SetCol(j + 1, tmp1);

        // Applying Givens Rotation to H col
        for (int i = 0; i <= j - 1; i++)
        {
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

        if (error <= tol)
            break;
    }

    // solve the triangular equation
    CSMatrix H_n(j + 1, j + 1);
    H.GetPartSlice(0, 0, j + 1, j + 1, H_n);

    CSMatrix Y(j + 1, 1);

    for (int i = 0; i <= j; i++)
    {
        Y(i, 0) = beta(i);
    }

    // Y.Print();
    Y.TRMM(H_n, 1.0, "L", "U");
    // Y.Print();

    // get solution x
    CSMatrix Q_n(n, j + 1);
    Q.GetPartSlice(0, 0, n, j + 1, Q_n);

    // Q.Print();
    CSVector y(j + 1);
    Y.GetColCopy(0, y);

    CSVector x_n(n);
    Q_n.GEMV(y, x_n);

    x += x_n;

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

        CSVector b(nv, b_data.data());

        CSVector x0(nv);

        x0 *= 0;

        double tol = 1e-12;
        int max_it = 500;

        auto sol_gmres = gmres(A, b, x0, tol, max_it);

        CSVector res(nv);
        A.SPMV(sol_gmres, res);

        res.Print();
        std::cout << res.Norm2() << std::endl;
        // sol_gmres.Print();

        delete xadj;
        delete adj;
        delete ew;
    }
    ChipSum::Common::Finalize();
}