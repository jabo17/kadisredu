//
// Created by jannickb on 6/17/24.
//
#include "kadisredu/data_structures/fast_set.h"

#include <kassert/kassert.hpp>

#include <numeric>

bool kadisredu::fast_set::contains(std::size_t i) const { return counter == set[i]; }

bool kadisredu::fast_set::add(std::size_t i) {
  KASSERT(i < set.size());
  if (!contains(i)) {
    set[i] = counter;
    return true;  // added
  }
  return false;  // already contained
}

void kadisredu::fast_set::clear() {
  if (counter == std::numeric_limits<int>::max()) {
    for (int &i : set) {
      i = 0;
      counter = 1;
    }
  } else {
    counter++;
  }
}
void kadisredu::fast_set::remove(std::size_t i) {
  KASSERT(i < set.size());
  set[i] = (counter - 1);
}
