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
    beta.HostToDevice();

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
            H.HostToDevice();

            Q.GetColCopy(i, tmp1);
            tmp1 *= (-1.0 * H(i, j));

            Q.GetColCopy(j + 1, tmp2);
            tmp2 += tmp1;

            Q.SetCol(j + 1, tmp2);
        }

        Q.GetColCopy(j + 1, tmp1);
        H(j + 1, j) = tmp1.Norm2();
        H.HostToDevice();

        tmp1 *= (1.0 / H(j + 1, j));
        Q.SetCol(j + 1, tmp1);

        // Applying Givens Rotation to H col
        for (int i = 0; i <= j - 1; i++)
        {
            CSFloat temp = cs(i) * H(i, j) + sn(i) * H(i + 1, j);
            H(i + 1, j) = -sn(i) * H(i, j) + cs(i) * H(i + 1, j);
            H(i, j) = temp;
        }
        H.HostToDevice();

        cs(j) = H(j, j) / sqrt(H(j, j) * H(j, j) + H(j + 1, j) * H(j + 1, j));
        cs.HostToDevice();

        sn(j) = H(j + 1, j) / sqrt(H(j, j) * H(j, j) + H(j + 1, j) * H(j + 1, j));
        sn.HostToDevice();

        H(j, j) = cs(j) * H(j, j) + sn(j) * H(j + 1, j);
        H(j + 1, j) = 0.0;
        H.HostToDevice();

        // update the residual vector
        beta(j + 1) = -sn(j) * beta(j);
        beta(j) = cs(j) * beta(j);
        beta.HostToDevice();

        error = abs(beta(j + 1)) / b.Norm2();
        printf("+++++++++++++++++++++++++++++++++++++++++++\n");
        printf("step # %d\n", j + 1);
        printf("residual : %.7f\n", error);
        
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
    Y.HostToDevice();

    Y.TRMM(H_n, 1.0, "L", "U");
    CSVector y(j + 1);
    Y.GetColCopy(0, y);

    CSMatrix Q_n(n, j + 1);
    Q.GetPartSlice(0, 0, n, j + 1, Q_n);

    // get solution x
    CSVector x_n(n);
    Q_n.GEMV(y, x_n);
    x += x_n;
    printf("+++++++++++++++++++++++++++++++++++++++++++\n");
    return x;
}

} // End namespace Solver
} // End namespace ChipSum