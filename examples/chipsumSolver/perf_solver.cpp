/* * * * * * * * * * * * * * * * * * * * *
 *   File:     cg-2.cpp
 *   Author:   Yaojie Yu
 *   group:    CDCS-HPC
 *   Time:     2022-06-12
 * * * * * * * * * * * * * * * * * * * * * */

#include <cg.hpp>
#include <bicg.hpp>
#include <bicgstab.hpp>
#include <gmres.hpp>
#include <iostream>

using namespace std;

#include <KokkosKernels_IOUtils.hpp>

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

        CSVector b(nv, b_data.data());

        CSVector x0(nv);

        x0 *= 0;

        x0 *= 0;

        double tol = 1e-5;
        int max_it = 100;

        Kokkos::Timer timer;
        // auto sol_bicg = ChipSum::Solver::cg(A, b, x0, tol, max_it);
        // auto sol_bicg = ChipSum::Solver::bicg(A, b, x0, tol, max_it);
        // auto sol_bicg = ChipSum::Solver::bicgstab(A, b, x0, tol, max_it);
        auto sol_bicg = ChipSum::Solver::gmres(A, b, x0, tol, max_it);
        double time = timer.seconds();

        cout << time << endl;
        delete xadj;
        delete adj;
        delete ew;
    }
    ChipSum::Common::Finalize();
}