#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_TYPES_H
#define TINYMESH_TYPES_H

#include <map>
#include <functional>

//! Pair of indices (i.e., unsigned int)
using IndexPair = std::pair<uint32_t, uint32_t>;

namespace std {
    //! Hash type for IndexPair
    template <>
    struct hash<IndexPair> {
        std::size_t operator()(const IndexPair &k) const {
            return std::get<0>(k) ^ std::get<1>(k);
        }
    };
}

#endif  // TINYMESH_TYPES_H
