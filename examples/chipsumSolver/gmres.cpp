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
    if (bnrm2 == 0.0) bnrm2 = 1.0;

    double error = r.Norm2() / bnrm2;

    if (error < tol) return x;

    Vector sn(m);
    Vector cs(m);
    Vector beta(m+1);

    beta(0) = r.Norm2();

    Matrix H(m+1, m);
    Matrix Q(n, m+1);

    int j = 0;
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

        auto sol_gmres= gmres(A, b, x0, tol, max_it);

        sol_gmres.Print();

        delete xadj;
        delete adj;
        delete ew;
    }
    ChipSum::Common::Finalize();
}