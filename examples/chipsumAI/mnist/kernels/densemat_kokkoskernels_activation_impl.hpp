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



namespace  ChipSum{
namespace  Numeric{
namespace  Impl {
namespace  DenseMat {

template <typename ValueType>
struct ReduceMaxFunctor{
    
    ValueType *_data;

    CHIPSUM_SPECIAL_INLINE
    ReduceMaxFunctor(ValueType *data) : _data(data) {}// printf("%f** ", 1.0);

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

        //_data(i,j) = _flagLog ？ ((_data(i,j) - row_max) - log(row_sum)) ： (exp(_data(i,j) - row_max) / row_sum);
    }
};


/* template <typename ValueType>
struct perRowMax_functor{
    
    Kokkos::View<ValueType **> _data;
    Kokkos::View<ValueType *> _out;

    HIPSUM_SPECIAL_INLINE
    perRowMax_functor(Kokkos::DualView<ValueType **> &data, Kokkos::View<ValueType *> &out) : _data(data.d_view), _out(out) {}

    CHIPSUM_SPECIAL_INLINE void operator() (const int64_t i, const int64_t j) const{
        _out(i) = max(_out(i), _data(i,j));
    }

};

template <typename ValueType>
struct perRowSum_functor{
    
    Kokkos::View<ValueType **> _data;
    Kokkos::View<ValueType *> _out;
    Kokkos::View<ValueType *> _row_max;

    HIPSUM_SPECIAL_INLINE
    perRowSum_functor(Kokkos::DualView<ValueType **> &data, Kokkos::View<ValueType *> &out, Kokkos::View<ValueType *> &row_max) : 
                        _data(data.d_view), _out(out), _row_max(row_max) {}

    CHIPSUM_SPECIAL_INLINE void operator() (const int64_t i, const int64_t j) const{
        _out(i) += exp(_data(i,j) - _row_max(i));
    }

};

template <typename ValueType>
struct Softmax_functor{
    
    Kokkos::View<ValueType **> _data;
    Kokkos::View<ValueType *> _out;

    Softmax_functor(Kokkos::DualView<ValueType **> &data, Kokkos::View<ValueType *> &out) : _data(data.d_view), _out(out) {}

    CHIPSUM_SPECIAL_INLINE void operator() (const int64_t i, const int64_t j) const{
        _out(i) += _data(i,j);
    }

}; */


template <typename ValueType>
CHIPSUM_FUNCTION_INLINE void
softmax(Kokkos::DualView<ValueType **> &A, bool flagLog = false) {
     const std::size_t M = A.extent(0);// row
     const std::size_t N = A.extent(1);// col
     
     typedef Kokkos::TeamPolicy<> team_policy;
     
     // flagLog控制是否进行logsoftmax计算
     SingleRowFunctor<ValueType> functor(A, N, flagLog);
     team_policy policy( M, N );
     Kokkos::parallel_for( "softmax", policy, functor);

     /* using mdrange_policy = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;

     Kokkos::View<ValueType *> row_max("pre row max val", M);
     perRowMax_functor<ValueType> functor_max(A, row_max);
     Kokkos::parallel_reduce( "get pre row max", mdrange_policy({0,0}, {M,N}), functor_max);

     Kokkos::View<ValueType *> row_sum("pre row sum val", M);
     perRowSum_functor<ValueType> functor_sum(A, row_sum, row_max);
     Kokkos::parallel_reduce( "get pre row sum", mdrange_policy({0,0}, {M,N}), functor_sum);

     Softmax_functor<ValueType> functor(A, row_sum, row_max);
     Kokkos::parallel_for( "softmax", mdrange_policy({0,0}, {M,N}), functor); */
}

}
}
}
}
#endif // DENSEMAT_KOKKOSKERNELS_ACTIVATION_IMPL_HPP
