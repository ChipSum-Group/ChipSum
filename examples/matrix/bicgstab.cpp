/* * * * * * * * * * * * * * * * * * * * *
 *   File:     bicg.cpp
 *   Author:   Yaojie Yu
 *   group:    CDCS-HPC
 *   Time:     2021-12-14
 * * * * * * * * * * * * * * * * * * * * * */

#include "../ChipSum.hpp"
#include "../chipsum/chipsum_macro.h"

Vector bicgstab(CSR &A, Vector &b, Vector &x, double tol, int max_it)
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

    double bnrm2 = b.Norm2();
    if (bnrm2 == 0.0)
        bnrm2 = 1.0;

    Vector r(x.GetSize());
    A.SPMV(x, r);
    b.AXPBY(r, 1.0, -1.0); // r = b - A*x

    double error = r.Norm2() / bnrm2;

    if (error < tol)
        return x;

    Vector r_tld(x.GetSize());
    r_tld.DeepCopy(r);

    Vector p(x.GetSize()), p_hat(x.GetSize());
    Vector v(x.GetSize()), t(x.GetSize());
    Vector s(x.GetSize()), s_hat(x.GetSize());

    double omega = 1.0;
    double alpha, beta, rho, rho_1;

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

        omega = (t.Dot(s)) / t.Dot(t);
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
        Vector r_temp(b.GetSize());
        A.SPMV(x, r_temp);          /* r_temp = A*x */
        b.AXPBY(r_temp, 1.0, -1.0); /* r_temp = b-r_temp */
        printf("%.20f\n", r_temp.Norm2());
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

        Vector b(nv, b_data.data());

        Vector x0(nv);

        x0 *= 0;

        double tol = 1e-12;
        int max_it = 500;

        auto sol_bicgstab = bicgstab(A, b, x0, tol, max_it);

        delete xadj;
        delete adj;
        delete ew;
    }
    ChipSum::Common::Finalize();
}
