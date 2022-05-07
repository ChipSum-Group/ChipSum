#include "../ChipSum.hpp"

int main(int argc, char *argv[]) {
    
    ChipSum::Common::Init(argc, argv);
    {   int N = 5;
        CSFloat v[5] = {1.0, 2.0, -6.1, 4.0, 5.0};
        CSVector a(N,v);

        ///norminf host输出接口
        CSFloat value = a.NormInf();
        std::cout<<"NormInf: "<<value<<std::endl;

        ///norminf device输出接口
        //目前unsigned才能通过Kokkos类型检测
        ChipSum::Numeric::Scalar<unsigned,ChipSum::Backend::DefaultBackend> scal(0);
        scal.Print();
        a.NormInf(scal);
        std::cout<<"NormInf: ";
        scal.Print();
    }
    ChipSum::Common::Finalize();
}
