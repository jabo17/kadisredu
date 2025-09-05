//
// Created by jannickb on 11/27/24.
//

#pragma once

#include "kadisredu/algorithms/mwis/greedy.h"
#include "kadisredu/algorithms/mwis/greedy_maximize.h"
#include "kadisredu/algorithms/mwis/mwis_def.h"

#include <kadisredu/algorithms/mwis/async_greedy.h>
#include <kadisredu/algorithms/mwis/async_greedy_maximize.h>
#include <kadisredu/algorithms/mwis/reducer.h>
#include <kadisredu/algorithms/mwis/sync_greedy.h>
#include <kadisredu/algorithms/mwis/sync_greedy_maximize.h>
#include <memory>

namespace kadisredu::mwis {

template <typename... Args>
std::unique_ptr<greedy> create_greedy_algo(const greedy_solver_context &ctx, Args &&...args) {
  using namespace kadisredu::mwis;
  if (ctx.greedy_algo == GreedyAlgo::SYNC) {
    return std::make_unique<sync::sync_greedy>(std::forward<Args>(args)...);
  } else {
    return std::make_unique<async::async_greedy>(std::forward<Args>(args)...,
                                                 ctx.message_queue_local_threshold);
  }
}

template <typename... Args>
std::unique_ptr<greedy_maximize> create_greedy_maximizer(ReducerType choice,
                                                         const Reducer::Config &cfg,
                                                         Args &&...args) {
  using namespace kadisredu::mwis;
  if (choice == ReducerType::SYNC) {
    return std::make_unique<sync::sync_greedy_maximize>(std::forward<Args>(args)...);
  } else {
    return std::make_unique<async::async_greedy_maximize>(std::forward<Args>(args)...,
                                                          cfg.message_queue_local_threshold);
  }
}

}  // namespace kadisredu::mwis
