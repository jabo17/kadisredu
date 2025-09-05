//
// Created by jannickb on 7/18/24.
//

#pragma once

#include <random>

namespace kadisredu {

class random_functions {
 public:
  static std::mt19937 mt;
  /**
   * Permutates a vector with of at least four elements fast (from KaHIP)
   * @tparam T
   * @param vec
   */
  template <typename T>
  static void fast_permutate(std::vector<T>& vec) {
    if(vec.size() <= 1) {
      return;
    }
    int distance = 20;
    std::uniform_int_distribution<int> A(0, distance);
    if (vec.size() < 4) {
      unsigned int posA = (A(random_functions::mt)) % vec.size();
      unsigned int posB = (A(random_functions::mt)) % vec.size();
      if(posA==posB) {
        posB = (posB + 1) % vec.size();
      }
      std::swap(vec[posA], vec[posB]);
      return;
    }
    unsigned int size = vec.size() - 4;
    for (unsigned int i = 0; i < size; i++) {
      unsigned int posA = i;
      unsigned int posB = (posA + A(random_functions::mt)) % size;
      std::swap(vec[posA], vec[posB]);
      std::swap(vec[posA + 1], vec[posB + 1]);
      std::swap(vec[posA + 2], vec[posB + 2]);
      std::swap(vec[posA + 3], vec[posB + 3]);
    }
  }

  static void reseed(int seed) { random_functions::mt.seed(seed); }

};
}  // namespace kadisredu