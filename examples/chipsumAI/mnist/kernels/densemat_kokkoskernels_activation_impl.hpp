///
/// \file     densemat_kokkoskernels_activation_impl.hpp
/// \author   Xiang Yukai
/// \group    CDCS-HPC
/// \date     
/// \brief    %stuff%
///

#ifndef __CHIPSUM_DENSEMAT_KOKKOSKERNELS_ACTIVATION_IMPL_HPP__
#define __CHIPSUM_DENSEMAT_KOKKOSKERNELS_ACTIVATION_IMPL_HPP__

#include <Kokkos_DualView.hpp>

#include "../chipsum/chipsum_macro.h"


template <typename ValueType>
struct ReduceMaxFunctor{
    
    ValueType *_data;

    CHIPSUM_SPECIAL_INLINE
    ReduceMaxFunctor(ValueType *data) : _data(data) {}

    CHIPSUM_SPECIAL_INLINE 
    void operator() (const int64_t j, ValueType & thisRowMax) const{
        thisRowMax = max(thisRowMax, _data[j]);
    }
};

template <typename ValueType>
struct ReduceSumFunctor{
    
    ValueType* _data;
    ValueType _row_max;

    CHIPSUM_SPECIAL_INLINE
    ReduceSumFunctor(ValueType *data, ValueType row_max) : 
                        _data(data), _row_max(row_max) {}

    CHIPSUM_SPECIAL_INLINE 
    void operator() (const int64_t j, ValueType & thisRowSum) const{
        // printf("%f*max* ", _row_max); // check reduce最大值是否正确
        thisRowSum += exp(_data[j] - _row_max);
    }
};

// softmax/logsoftmax 行处理functor
template <typename ValueType>
struct SingleRowFunctor{

    typedef Kokkos::TeamPolicy<>::member_type  member_type;
    Kokkos::View<ValueType **> _data;
    bool _flagLog;
    const int64_t _N;

    SingleRowFunctor(Kokkos::DualView<ValueType **> &data, const int64_t N, bool flagLog) : 
                        _data(data.d_view), _N(N), _flagLog(flagLog) {}

    CHIPSUM_SPECIAL_INLINE 
    void operator() (const member_type & teamMember) const{
        const int i = teamMember.league_rank();
        ValueType row_max = _data(i,0);
        ValueType row_sum = 0.0;

        ReduceMaxFunctor<ValueType> functor_max(&(_data.data()[i*_N]));// 第i行的首地址
        ReduceSumFunctor<ValueType> functor_sum(&(_data.data()[i*_N]), row_max);// 第i行的首地址, 及该行最大值
        
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange( teamMember, _N ), functor_max, row_max);
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange( teamMember, _N ), functor_sum, row_sum);

        const int j = teamMember.team_rank();
        if(_flagLog) _data(i,j) =(_data(i,j) - row_max) - log(row_sum);// logSoftmax
        else _data(i,j) = exp(_data(i,j) - row_max) / row_sum;// Softmax

    }
};



template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
softmax(ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &input, bool flagLog = false) {
     auto A = input.GetData();
     const std::size_t M = A.extent(0);// row
     const std::size_t N = A.extent(1);// col
     
     typedef Kokkos::TeamPolicy<> team_policy;
     
     // flagLog控制是否进行logsoftmax计算
     SingleRowFunctor<ValueType> functor(A, N, flagLog);
     team_policy policy( M, N );
     Kokkos::parallel_for( "softmax", policy, functor);

}

template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void argmax(ChipSum::Numeric::DenseMatrix<ValueType, ChipSum::Backend::DefaultBackend> &input) {
    auto A = input.GetData();
    
    ::std::size_t M = A.extent(0);
    ::std::size_t N = A.extent(1);
    // auto h_A = Kokkos::create_mirror_view(A);

    Kokkos::deep_copy(A.h_view, A.d_view);

    // cout << "prediction is" << ":" << endl;

    ValueType max_position = 0;
    ValueType max_val = A.h_view(0, 0);

    for (std::size_t i = 0; i < N; ++i) {
        if(A.h_view(0, i)>=max_val){
            max_val = A.h_view(0, i);
            max_position = i;
        }
    }
    ::std::cout << "*****prediction is*****" << " : " << max_position << ::std::endl;
    ::std::cout << ::std::endl;
}


#endif // DENSEMAT_KOKKOSKERNELS_ACTIVATION_IMPL_HPP
