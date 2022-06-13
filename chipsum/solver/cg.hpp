/* * * * * * * * * * * * * * * * * * * * *
 *   File:     cg.cpp
 *   Author:   Yaojie Yu
 *   group:    CDCS-HPC
 *   Time:     2022-06-07
 * * * * * * * * * * * * * * * * * * * * * */

#include "../../ChipSum.hpp"
#include "../chipsum_macro.h"

namespace ChipSum {
namespace Solver {
void cg(CSR &A, CSVector &b, CSVector &x, CSFloat tol, int max_it)
{
    // Conjugate Gradient Method without preconditioning.
    //
    // input   A        REAL matrix
    //         x        REAL initial guess vector
    //         b        REAL right hand side vector
    //         tol      REAL error tolerance
    //         max_it   INTEGER maximum number of iterations
    //
    // output  x        REAL solution vector

    CSVector r(x.GetSize());
    A.SPMV(x, r);
    b.AXPBY(r, 1.0, -1.0); // r = b - A*x

    CSVector p(b.GetSize());
    p.DeepCopy(r);

    CSVector Ap(b.GetSize());

    CSFloat alpha = 0, beta = 0.0, rsnew = 0;
    CSFloat rsold = r.Dot(r);
    CSFloat error;

    for (int i = 0; i < max_it; i++)
    {

        A.SPMV(p, Ap);

        alpha = rsold / (p.Dot(Ap));

        p.AXPBY(x, alpha, 1.0);

        Ap.AXPBY(r, -alpha, 1.0);

        rsnew = r.Dot(r);
        
        error = sqrt(rsnew);
        
        if (error < tol){
            printf("+++++++++++++++++++++++++++++++++++++++++++\n");
            printf("step # %d\n", i + 1);
            printf("residual : %.7f\n", error);
            break;
        }

        beta = rsnew / rsold;
        r.AXPBY(p, 1.0, beta);

        rsold = rsnew;
    }
    printf("+++++++++++++++++++++++++++++++++++++++++++\n");
}

} // End namespace Solver
} // End namespace ChipSum