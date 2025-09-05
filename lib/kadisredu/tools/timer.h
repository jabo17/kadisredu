//
// Created by jannickb on 6/26/24.
// Based on timer.h from KaMinPar
//

#pragma once

#include <chrono>
#include <string_view>

namespace kadisredu {

#define KADISREDU_SCOPED_TIMER2(handle) \
  timer t_##__LINE__;                   \
  auto sc_##__LINE__ = scoped_timer(&t_##__LINE__, handle);

#define KADISREDU_SCOPED_TIMER1(timer_ptr) KADISREDU_SCOPED_TIMER2(timer_ptr, []() {})

#define KADISREDU_DEFINE_METHOD_TIMER(ident) method_timer ident;

#define KADISREDU_MEASURE_TIME_METHOD(ident) auto s_t##__LINE__ = ident.measure_method_call();

class timer {
 public:
  using time_point = std::chrono::time_point<std::chrono::high_resolution_clock>;
  using duration = std::chrono::high_resolution_clock::duration;

  timer() = default;

  static time_point now() { return std::chrono::high_resolution_clock::now(); }
  void restart() { m_start = timer::now(); }
  double elapsed_seconds() { return seconds(timer::now() - m_start); }
  double final_elapsed_seconds() { return seconds(m_total_elapsed); }
  void stop() { m_total_elapsed = timer::now() - m_start; }
  static double seconds(duration dur) {
    return static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(dur).count()) /
           1000.0;
  }

 private:
  time_point m_start{};
  duration m_total_elapsed{};
};

template <typename StopHandle>
class scoped_timer {
 public:
  explicit scoped_timer(timer* t, StopHandle on_stop) : m_t(t), m_on_stop(on_stop) {
    m_t->restart();
  }

  inline ~scoped_timer() {
    m_t->stop();
    m_on_stop();
  }

 private:
  timer* m_t;
  StopHandle m_on_stop;
};

class method_timer {
 public:
  method_timer() = default;

  auto measure_method_call() {
    return scoped_timer(&m_method_timer,
                        [this]() { m_seconds += m_method_timer.final_elapsed_seconds(); });
  }

  [[nodiscard]] double total_seconds() const { return m_seconds; }

 protected:
  double m_seconds = 0;
  timer m_method_timer;
};

}  // namespace kadisredu
