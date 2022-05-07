#include "../ChipSum.hpp"

int main(int argc, char *argv[]) {
    
    ChipSum::Common::Init(argc, argv);
    {   int N=10;
        CSFloat *v = static_cast<CSFloat *>(std::malloc(N * sizeof(CSFloat)));
        for (int i = 0; i < N; ++i) {
            v[i] = 1.0;
        }
        
        CSVector a(N,v);
        // CSFloat alfa=2.0;
        // a.Fill(alfa);
        a.Fill(2.0);
        std::cout<<"fill vector";
        a.Print();     
        std::free(v);
    }
    ChipSum::Common::Finalize();
}
