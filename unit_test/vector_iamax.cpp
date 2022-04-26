#include "../ChipSum.hpp"

int main(int argc, char *argv[]) {
    
    ChipSum::Common::Init(argc, argv);
    {   int N = 5;
        CSFloat v[5] = {1.0, 2.0, -6.0, 6.0, 4.0};
        CSVector a(N,v);

        ///IAMAX host输出接口
        size_t loc = a.IAMAX();
        std::cout<<"iamax1 location: "<<loc<<std::endl;

        ///IAMAX device输出接口
        //目前unsigned才能通过Kokkos类型检测
        ChipSum::Numeric::Scalar<unsigned,ChipSum::Backend::DefaultBackend> scal(0);
        std::cout<<"iamax2 location origin: ";
        scal.Print();
        a.IAMAX(scal);
        std::cout<<"iamax2 location: ";
        scal.Print();
    }
    ChipSum::Common::Finalize();
}
