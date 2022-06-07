/* * * * * * * * * * * * * * * * * * * * *
 *   File:     bicg.cpp
 *   Author:   Yaojie Yu
 *   group:    CDCS-HPC
 *   Time:     2022-06-07
 * * * * * * * * * * * * * * * * * * * * * */

#include "../../ChipSum.hpp"
#include "../chipsum_macro.h"

namespace ChipSum {
namespace Solver {
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

} // End namespace Solver
} // End namespace ChipSum