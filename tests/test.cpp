///
/// \file     test.cpp
/// \author   Riiiichman-Li
/// \group    CDCS-HPC
/// \date     2021-12-15
/// \brief    %stuff%
///

#include <iostream>
using namespace std;

#include <KokkosKernels_IOUtils.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosKernels_default_types.hpp>



#include <type_traits>
#include <vector>




#include "../ChipSum.hpp"
#include "../chipsum/chipsum_macro.h"



using Scal  = default_scalar;
using Ordinal = default_lno_t;
using Offset  = default_size_type;
using Layout  = default_layout;


CHIPSUM_FUNCTION_INLINE void run_cg(CSR& A,CSVector& b,CSVector& x,double tol=10e-12,int max_it=200)
{

    //    x0 = np.zeros(len(b))
    //    r0 = b-np.dot(A,x0)
    //    p0 = r0
    CSVector x0(x.GetSize());

    CSVector r0(x.GetSize());

    A.SPMV(x0,r0);

    b.AXPBY(r0,1.,-1.);

    CSVector p0(x.GetSize());

    p0.DeepCopy(r0);

    CSVector Ap(x.GetSize());

    CSVector err(x.GetSize());


    for(int i=0;i<max_it;++i){
        //        alpha = np.dot(r0.T,r0)/np.dot(p0.T,np.dot(A,p0.T))
        double alpha=r0.Dot(r0);

        A.SPMV(p0,Ap);
        alpha /= p0.Dot(Ap);

        //        x1 = x0+alpha*p0
        x.DeepCopy(x0);
        p0.AXPBY(x,alpha);

        //        r1 = r0-alpha*np.dot(A,p0)
        CSVector r1(r0.GetSize());
        r1.DeepCopy(r0);
        Ap.AXPBY(r1,-alpha,1);

        //        beta = np.dot(r1.T,r1)/np.dot(r0.T,r0)
        double beta = r1.Dot(r1);
        beta /= r0.Dot(r0);

        //        p0 = r1+beta*p0
        r1.AXPBY(p0,1.,beta);

        //        x0 = x1;r0 = r1
        x0.DeepCopy(x);
        r0.DeepCopy(r1);

        A.SPMV(x0,err);
        b.AXPBY(err,1,-1);

        printf("%.28f\n",err.Norm2());

        if(err.Norm2()<tol) {


            return;
        }


    }


}
typedef ChipSum::Numeric::DenseMatrix<CSFloat,ChipSum::Backend::Serial>
SerialMatrix;
inline void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn)
{
   double temp = cs * dx + sn * dy;
   dy = -sn * dx + cs * dy;
   dx = temp;
}


int main(int argc, char *argv[]) {

    const char* filename_A = argv[1];
    const char* filename_b = argv[2];
    if(filename_A == nullptr) filename_A = "../../data/A.mtx";

    ChipSum::Common::Init(argc, argv);
    {

        CSInt nv = 0, ne = 0;
        CSInt *xadj, *adj;
        double *ew;

        KokkosKernels::Impl::read_matrix<CSInt,CSInt, double> (&nv, &ne, &xadj, &adj, &ew, filename_A);

        CSR A(nv,nv,ne,xadj,adj,ew);

        A.SavePatternFig("../A.PNG");


        vector<double> b_data;
        double temp;

        ifstream IN(filename_b);

        for(int i=0;i<nv;++i){
            temp = 1.0;
            b_data.push_back(temp);
        }


        IN.close();

        CSVector b(nv,b_data.data());

        CSVector x(nv);



        run_cg(A,b,x,10e-5,1477);

        delete xadj;
        delete adj;
        delete ew;
    }
    ChipSum::Common::Finalize();
}