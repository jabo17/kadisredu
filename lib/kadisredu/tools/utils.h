//
// Created by jannickb on 7/4/24.
//

#pragma once

#include "kadisredu/data_structures/fast_set.h"
#include "kadisredu/data_structures/sized_vector.h"

#include <span>

namespace kadisredu::utils {

template <typename InputIt>
constexpr inline std::size_t intersection(const fast_set& set, InputIt first, InputIt last) {
  std::size_t result = 0;

  for (; first != last; ++first) {
    if (set.contains(*first)) {
      ++result;
    }
  }

  return result;
}

template <typename InputIt>
constexpr inline bool disjoint(const fast_set& set, InputIt first, InputIt last) {
  for (; first != last; ++first) {
    if (set.contains(*first)) {
      return false;
    }
  }
  return true;
}

template <typename InputIt>
void inline add_to_set(fast_set& set, InputIt first, InputIt last) {
  for (; first != last; ++first) {
    set.add(*first);
  }
}

/**
    @brief partitions a sequence in two sequences by reordering them according to a predicate.
    The first partition contains all elements not matching the predicate and the second contains
   those that match the predicate.
    @return an itarator to the first element matching the predicate; otherwise @param last
*/
template <typename Container, typename UnaryPred>
Container::iterator partition(Container& a, UnaryPred&& p) {
  using std::swap;

  if (a.size() == 0) {
    return a.end();
  }
  long long i = 0;
  long long j = a.size() - 1;
  do {
    while (i < a.size() && !p(a[i])) ++i;
    while (j >= 0 && p(a[j])) --j;
    if (i < j) {
      swap(a[i], a[j]);
    }
  } while (i <= j);
  KASSERT(a.begin() + i == a.end() || p(a[i]));
  KASSERT(i == 0 || !p(a[i - 1]));
  return a.begin() + i;
}

}  // namespace kadisredu::utils
