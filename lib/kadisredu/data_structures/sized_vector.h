//
// Created by jannickb on 6/25/24.
//

#pragma once

#include <kassert/kassert.hpp>
#include <memory>
#include <span>
#include <vector>
namespace kadisredu {

template <typename T>
class sized_vector {
 public:
  using size_t = std::size_t;
  using iterator = std::vector<T>::iterator;
  using const_iterator = std::vector<T>::const_iterator;
  using value_type = T;

  explicit sized_vector(size_t capacity) : m_data(capacity), counter(0) {}

  sized_vector() : m_data(0), counter(0) {}

  void push_back(T&& element) {
    KASSERT(counter < m_data.size());
    m_data[counter++] = std::move(element);
  }

  void push_back(const T& element) {
    KASSERT(counter < m_data.size());
    m_data[counter++] = element;
  }

  void pop_back() {
    KASSERT(counter > 0);
    --counter;
  }

  void pop_back(std::size_t elements) {
    KASSERT(counter >= elements);
    counter -= elements;
  }

  [[nodiscard]] T& back() {
    KASSERT(counter > 0);
    return m_data[counter - 1];
  }

  [[nodiscard]] T& front() {
    KASSERT(counter > 0);
    return m_data[0];
  }
  void clear() { counter = 0; }

  [[nodiscard]] size_t capacity() { return m_data.size(); }

  [[nodiscard]] size_t size() { return counter; }

  [[nodiscard]] size_t size() const { return counter; }

  void resize(std::size_t new_size) {
    KASSERT(new_size <= capacity());
    counter = new_size;
  }

  // explicit operator std::span<T>() { return std::span<T>(data(), size()); }
  // explicit operator std::span<const T>() const { return std::span<const T>(data(), size()); }

  T* data() { return &m_data[0]; }
  const T* data() const { return &m_data[0]; }

  bool empty() { return size() == 0; }

  [[nodiscard]] constexpr iterator begin() { return m_data.begin(); }

  [[nodiscard]] constexpr const_iterator begin() const { return m_data.begin(); }

  [[nodiscard]] constexpr iterator end() { return m_data.begin() + counter; }

  [[nodiscard]] constexpr const_iterator end() const { return m_data.begin() + counter; }

  [[nodiscard]] T& operator[](size_t pos) {
    KASSERT(pos < counter);
    return m_data[pos];
  }

  [[nodiscard]] const T& operator[](size_t pos) const {
    KASSERT(pos < counter);
    return m_data[pos];
  }

 private:
  std::vector<T> m_data;
  size_t counter;
};

}  // namespace kadisredu
