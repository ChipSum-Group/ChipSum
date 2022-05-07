#include "../ChipSum.hpp"

int main(int argc, char *argv[]) {
    
    ChipSum::Common::Init(argc, argv);
    {   int N = 5;
        CSFloat v[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
        CSVector a(N,v);
        CSFloat sum = a.Sum();
        std::cout<<"vector sum : "<<sum<<std::endl;

    }
    ChipSum::Common::Finalize();
}
