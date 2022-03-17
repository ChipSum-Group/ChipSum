
/* * * * * * * * * * * * * * * * * * * * *
 *   File:     test3.cpp
 *   Author:   Yaojie Yu
 *   group:    HPC group@NSCC-CD
 *   Time:     2022-03-18
 * * * * * * * * * * * * * * * * * * * * * */

#include "../ChipSum.hpp"
#include "../chipsum/chipsum_macro.h"

#define K 2
#define B 2
#define N 5
#define M 3

int main(int argc, char *argv[])
{
    ChipSum::Common::Init(argc, argv);
    {
        double *m1 = static_cast<double *>(std::malloc(5 * 5 * sizeof(double)));
        for (int i = 0; i < 5 * 5; ++i)
        {
            m1[i] = double(1);
        }

        Matrix M1(5, 5, m1);
        M1.Print();

        double *m2 = static_cast<double *>(std::malloc(5 * sizeof(double)));
        for (int i = 0; i < 5; ++i)
        {
            m2[i] = double(i);
        }

        Vector M2(5, m2);
        M2.Print();

        M1.SetCol(1, M2);
        M1.Print();
        M1.SetRow(1, M2);
        M1.Print();
    }
    ChipSum::Common::Finalize();
}
