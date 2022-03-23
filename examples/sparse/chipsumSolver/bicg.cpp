/* * * * * * * * * * * * * * * * * * * * *
 *   File:     bicg.cpp
 *   Author:   Yaojie Yu
 *   group:    CDCS-HPC
 *   Time:     2021-12-14
 * * * * * * * * * * * * * * * * * * * * * */

#include "../ChipSum.hpp"
#include "../chipsum/chipsum_macro.h"

CSVector bicg(CSR &A, CSVector &b, CSVector &x, double tol, int max_it)
{
    // BiConjugate Gradient Method without preconditioning.
    //
    // input   A        REAL matrix
    //         x        REAL initial guess vector
    //         b        REAL right hand side vector
    //         tol      REAL error tolerance
    //         max_it   INTEGER maximum number of iterations
    //
    // output  x        REAL solution vector

    double bnrm2 = b.Norm2();
    if (bnrm2 == 0.0)
        bnrm2 = 1.0;

    CSVector r(x.GetSize());
    A.SPMV(x, r);
    b.AXPBY(r, 1.0, -1.0); // r = b - A*x

    double error = r.Norm2() / bnrm2;

    if (error < tol)
        return x;

    CSVector r_tld(x.GetSize());
    r_tld.DeepCopy(r);

    CSVector p(x.GetSize()), p_tld(x.GetSize());
    CSVector z(x.GetSize()), z_tld(x.GetSize());
    CSVector q(x.GetSize()), q_tld(x.GetSize());

    double alpha, beta, rho, rho_1;

    for (int i = 0; i < max_it; i++)
    {

        // z = M.inverse() *r; M is a precondtioner, here M is a identity matrix
        // z_tld = M'.inverse() *r_tld; M is a precondtioner

        z.DeepCopy(r);
        z_tld.DeepCopy(r_tld);

        rho = r.Dot(z_tld);

        if (rho == 0.0)
            break;

        if (i > 0)
        {
            beta = rho / rho_1;
            z.AXPBY(p, 1.0, beta);         // p = z + beta * p;
            z_tld.AXPBY(p_tld, 1.0, beta); // p_tld = z_tld + beta * p_tld;
        }
        else
        {
            p.DeepCopy(z);         // p = z
            p_tld.DeepCopy(z_tld); // p_tld = z_tld
        }

        A.SPMV(p, q);         // q = A*p
        A.SPMV(p_tld, q_tld); // q_tld = A'x, but here A is symmetry

        alpha = rho / (p_tld.Dot(q));

        p.AXPBY(x, alpha, 1.0);          // x = x + alpha * p;
        q.AXPBY(r, -alpha, 1.0);         // r = r - alpha * Ap
        q_tld.AXPBY(r_tld, -alpha, 1.0); // r_tld = r_tld - alpha * q_tld

        error = r.Norm2() / bnrm2;

        if (error <= tol)
            break;

        CSVector r_temp(b.GetSize());
        A.SPMV(x,r_temp); /* r_temp = A*x */
        b.AXPBY(r_temp,1.0,-1.0); /* r_temp = b-r_temp */
        printf("%.20f\n",r_temp.Norm2());

        rho_1 = rho;
    }

    return x;
}

int main(int argc, char *argv[])
{

    char *filename_A = argv[1];
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

        CSVector b(nv,b_data.data());

        CSVector x0(nv);

        x0 *= 0;

        double tol = 1e-12;
        int max_it = 500;

        auto sol_bicg = bicg(A, b, x0, tol, max_it);

        delete xadj;
        delete adj;
        delete ew;
    }
    ChipSum::Common::Finalize();
}
