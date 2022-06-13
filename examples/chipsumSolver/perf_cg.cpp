/* * * * * * * * * * * * * * * * * * * * *
 *   File:     cg-2.cpp
 *   Author:   Yaojie Yu
 *   group:    CDCS-HPC
 *   Time:     2022-06-12
 * * * * * * * * * * * * * * * * * * * * * */

#include <cg.hpp>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{

    ChipSum::Common::Init(argc, argv);
    {
        CSInt N = 20000;

        CSInt nnz     = 2 + 3*(N - 2) + 2;
        CSInt* ptrRaw = new CSInt[N + 1];
        CSInt* indRaw = new CSInt[ nnz ];
        CSFloat* valRaw  = new CSFloat[ nnz ];


    // Add rows one-at-a-time
    for (int i = 0; i < (N + 1); i++) {
        if (i==0) {
            ptrRaw[0] = 0;
            indRaw[0] = 0;   indRaw[1] = 1;
            valRaw[0] = 2.0; valRaw[1] = -1.0;
        }
        else if (i==N) {
            ptrRaw[N] = nnz;
        }
        else if (i==(N-1)) {
            ptrRaw[i] = 2 + 3*(i-1);
            indRaw[2 + 3*(i-1)] = i-1;  indRaw[2 + 3*(i-1) + 1] = i;
            valRaw[2 + 3*(i-1)] = -1.0; valRaw[2 + 3*(i-1) + 1] = 2.0;
        }
        else {
            ptrRaw[i] = 2 + 3*(i-1);
            indRaw[2 + 3*(i-1)] = i-1;  indRaw[2 + 3*(i-1) + 1] = i;   indRaw[2 + 3*(i-1) + 2] = i+1;
            valRaw[2 + 3*(i-1)] = -1.0; valRaw[2 + 3*(i-1) + 1] = 2.0; valRaw[2 + 3*(i-1) + 2] = -1.0;
        }
    }

         CSR A(N,N,nnz,ptrRaw,indRaw,valRaw);

        //  A.Print();

        // CSR A(nv, nv, ne, xadj, adj, ew);


        vector<double> b_data;
        double temp;

        // for (int i = 0; i < nv; ++i)
        for (int i = 0; i < N; ++i)
        {
            temp = 1.;
            b_data.push_back(temp);
        }

        // CSVector b(nv, b_data.data());
        // CSVector x0(nv);
        CSVector b(N, b_data.data());
        CSVector x0(N);
        

        x0 *= 0;

        double tol = 1e-5;
        int max_it = 100000;
        ChipSum::Solver::cg(A, b, x0, tol, max_it);

        Kokkos::Timer timer;
        for (int i=0; i<10; ++i)
            {x0 *= 0;
             ChipSum::Solver::cg(A, b, x0, tol, max_it);
            }
        double time = timer.seconds();

        cout << time << endl;
        delete ptrRaw;
        delete indRaw;
        delete valRaw;
    }
    ChipSum::Common::Finalize();
}