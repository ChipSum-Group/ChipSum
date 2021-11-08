#ifndef __CHIPSUM_DATA_TYPES_HPP__
#define __CHIPSUM_DATA_TYPES_HPP__

#include <cctype>

namespace ChipSum {
struct DataTypes{
    using float64 = double;
    using float32 = float;

    /// \brief AKA: int
    using int32 = __int32_t;
    /// \brief AKA: unsigned int
    using uint32 = __uint32_t;

};
}

#endif
