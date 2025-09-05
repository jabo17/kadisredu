//
// Created by Jannick Borowitz on 6/17/24.
// Based on fast_set.h from KaHIP
//
#pragma once

#include <vector>

namespace kadisredu {
class fast_set {
 public:
  explicit fast_set(std::size_t n) : set(n, 0), counter(1) {}

  [[nodiscard]] bool contains(std::size_t i) const;
  bool add(std::size_t i);
  void remove(std::size_t i);
  void clear();

  void resize(std::size_t n) {
    set.resize(n, 0);
    clear();
  }

 private:
  std::vector<int> set;
  int counter = 0;
};
};  // namespace kadisredu
