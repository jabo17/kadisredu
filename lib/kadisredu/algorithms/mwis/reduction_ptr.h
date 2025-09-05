//
// Created by jannickb on 10/21/24.
//

#pragma once

#include "kadisredu/definitions.h"
#include <memory>

namespace kadisredu {

template <typename Reduction>
class ReductionPtr {
 private:
  std::unique_ptr<Reduction> ptr_;

 public:
  ReductionPtr() : ptr_(nullptr) {}

  explicit ReductionPtr(Reduction* reduction) noexcept : ptr_(reduction) {}

  ReductionPtr(const ReductionPtr& other) : ptr_(std::move(other.ptr_->clone())) {}

  ReductionPtr(ReductionPtr&& other) noexcept : ptr_(std::move(other.ptr_)) {
    other.ptr_ = nullptr;
  }

  ReductionPtr& operator=(const ReductionPtr& other) {
    if (this == &other) {
      return *this;
    }
    ptr_ = std::move(other.ptr_->clone());
    return *this;
  }

  ReductionPtr& operator=(ReductionPtr&& other) noexcept {
    ptr_ = std::move(other.ptr_);
    return *this;
  }

  Reduction& operator*() { return *ptr_; }
  Reduction* operator->() const noexcept { return ptr_.get(); }
};

}