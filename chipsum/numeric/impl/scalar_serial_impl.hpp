/*
 * @Description: 标量scalar的串行实现
 * @Version: 2.0
 * @Autor: Li Kunyun
 * @Date: 2021-08-09 13:58:38
 * @LastEditors: Li Kunyun
 * @LastEditTime: 2021-08-12 10:39:52
 */


#ifndef __CHIPSUM_SCALAR_SERIAL_IMPL_HPP__
#define __CHIPSUM_SCALAR_SERIAL_IMPL_HPP__

#include <fstream>
#include <vector>

#include "../../backend/backend.hpp"
#include "../../chipsum_macro.h"
#include "../numeric_traits.hpp"

namespace ChipSum {
namespace Numeric {

template <typename ScalarType, typename SizeType, typename... Props>
struct Scalar_Traits<ScalarType, SizeType, ChipSum::Backend::Serial, Props...>
    : public Operator_Traits<ScalarType, SizeType, ChipSum::Backend::Serial> {
  using scalar_type = ScalarType;

  using device_scalar_value_type = ScalarType;
};

namespace Impl {
namespace Scalar {
template <typename ScalarType, typename SizeType, typename... Props>
CHIPSUM_FUNCTION_INLINE void Create(ScalarType &) {
  return;
}

template <typename ScalarType, typename SizeType, typename... Props>
CHIPSUM_FUNCTION_INLINE void Create(const ScalarType s, ScalarType &r) {
  r = s;
}

template <typename ScalarType, typename SizeType, typename... Props>
CHIPSUM_FUNCTION_INLINE void DeepCopy(const ScalarType s, ScalarType &r) {
  r = s;
}

template <typename ScalarType, typename SizeType, typename... Props>
CHIPSUM_FUNCTION_INLINE ScalarType GetItem(ScalarType s) {
  return s;
}

template <typename ScalarType, typename SizeType, typename... Props>
CHIPSUM_FUNCTION_INLINE void GetItem(const ScalarType s, ScalarType &r) {
  r = s;
}

template <typename ScalarType, typename SizeType, typename... Props>
CHIPSUM_FUNCTION_INLINE void Mult(const ScalarType s,
                                  const std::vector<ScalarType> &v,
                                  std::vector<ScalarType> &r) {
  assert(v.size() == r.size());
  for (std::size_t i = 0; i < r.size(); ++i) {
    r[i] = s * v[i];
  }
}
template <typename ScalarType, typename SizeType, typename... Props>
CHIPSUM_FUNCTION_INLINE void Print(const ScalarType s, std::ostream &out) {
  cout << s << endl;
}
} // namespace Scalar
} // namespace Impl

} // namespace Numeric
} // namespace ChipSum

#endif