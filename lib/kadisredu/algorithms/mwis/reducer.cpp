//
// Created by jannickb on 12/31/24.
//

#include "reducer.h"

#include "kadisredu/algorithms/mwis/weighted_reductions.h"

namespace kadisredu::mwis {
void Reducer::init_reduction_style(Reducer::ReduceStatus &status) {
  switch (status.cfg.style) {
    case ReductionStyle::PRELIMINARY:
      status.reductions = make_reduction_vector<neighborhood_removal, degree_one>(
          status.graph.num_nodes(), status.graph.num_ghosts());
      break;
    case ReductionStyle::NORMAL:
      status.reductions =
          make_reduction_vector<degree_one, neighborhood_removal, simplicial_weight_transfer,
                                v_shape, basic_single_edge, extended_single_edge,
                                neighborhood_folding, generalized_neighborhood_removal,
                                generalized_neighborhood_folding>(  //
              status.graph.num_nodes(), status.graph.num_ghosts());
      break;
    case ReductionStyle::REDUCE_AND_PEEL:
      status.reductions =
          make_reduction_vector<degree_one, neighborhood_removal, simplicial_weight_transfer,
                                v_shape, basic_single_edge, extended_single_edge,
                                neighborhood_folding, generalized_neighborhood_removal,
                                generalized_neighborhood_folding, peeler>(  //
              status.graph.num_nodes(), status.graph.num_ghosts());
      break;
    case ReductionStyle::FULL:
    default:
      status.reductions = make_reduction_vector<
          degree_one, neighborhood_removal, clique, simplicial_weight_transfer, v_shape, domination,
          basic_single_edge, extended_single_edge, neighborhood_folding,
          generalized_neighborhood_removal, generalized_neighborhood_folding, extended_domination>(  //
          status.graph.num_nodes(), status.graph.num_ghosts());
  }
  status.reduction_map.resize(MAX_REDUCTION_NUM + 1);
  for (std::size_t i = 0; i < status.reductions.size(); ++i) {
    status.reduction_map[static_cast<std::size_t>(status.reductions[i]->get_reduction_type())] = i;
  }
  status.reduction_measurements = std::vector<ReduceMeasurement>(
      status.reductions.size() + 1);  // one additional stat for received u
}

void Reducer::mark_node(NodeID node) {
  KASSERT(status_.graph.is_visible(node) && !status_.graph.is_ghost(node));
  for (auto &reduction : status_.reductions) {
    if (!reduction->marker.mark(node)) {
      KASSERT(status_.reductions.back()->marker.marked_next_set.contains(node));
      return;
    }
  }
}
void Reducer::mark_neigh(NodeID node) {
  for (auto &target : status_.graph[node]) {
    if (!status_.graph.is_ghost(target)) {
      mark_node(target);
    }
  }
}

std::optional<Reducer::KernelStatistic> Reducer::aggregate_reduce_statistic() {
  using namespace kamping;

  struct reduce_statistic_message {
    // overall (sum)
    GlobalNodeID kernel_nodes;
    GlobalNodeWeight offset;

    std::size_t allreduce;
    std::size_t alltoallv;
    std::size_t max_recv_border_updates;
    std::size_t sum_recv_border_updates;
    std::size_t sum_recv_redundant_border_updates;
    double max_recv_redundant_border_updates_rel;
  };

  reduce_statistic_message stat_msg = {
      .offset = status_.get_solution_weight(),
      .allreduce = status_.allreduce,
      .alltoallv = status_.alltoallv,
      .max_recv_border_updates = status_.recv_border_updates,
      .sum_recv_border_updates = status_.recv_border_updates,
      .sum_recv_redundant_border_updates = status_.recv_redundant_border_updates,
      .max_recv_redundant_border_updates_rel =
          static_cast<double>(status_.recv_redundant_border_updates) /
          static_cast<double>(status_.recv_border_updates)};

  auto agg_stat = comm_.reduce_single(
      send_buf(stat_msg),
      op(
          [](auto &s1, auto &s2) {
            return reduce_statistic_message{
                .offset = s1.offset + s2.offset,
                .allreduce = std::max(s1.allreduce, s2.allreduce),
                .alltoallv = std::max(s1.alltoallv, s2.alltoallv),
                .max_recv_border_updates =
                    std::max(s1.max_recv_border_updates, s2.max_recv_border_updates),
                .sum_recv_border_updates = s1.sum_recv_border_updates + s2.sum_recv_border_updates,
                .sum_recv_redundant_border_updates =
                    s1.sum_recv_redundant_border_updates + s2.sum_recv_redundant_border_updates,
                .max_recv_redundant_border_updates_rel =
                    std::max(s1.max_recv_redundant_border_updates_rel,
                             s2.max_recv_redundant_border_updates_rel)};
          },
          ops::commutative));

  auto kernel_graph_stat = graph_summary::agg_graph_summary_on_root(comm_, status_.graph);

  struct reduction_statistic_message {
    NodeID min_reduced_local_nodes;
    NodeID max_reduced_local_nodes;
    NodeID min_reduced_ghosts;
    NodeID max_reduced_ghosts;
  };

  std::vector<reduction_statistic_message> reduction_statistics(status_.reductions.size() + 1);
  for (std::size_t i = 0; i < status_.reduction_measurements.size(); ++i) {
    auto &stat = status_.reduction_measurements[i];

    reduction_statistics[i] = {.min_reduced_local_nodes = stat.reduced_local_nodes,
                               .max_reduced_local_nodes = stat.reduced_local_nodes,
                               .min_reduced_ghosts = stat.reduced_ghosts,
                               .max_reduced_ghosts = stat.reduced_ghosts};
  }

  auto agg_reduction_statistics = comm_.reduce(
      send_buf(reduction_statistics),
      op(
          [](auto &s1, auto &s2) {
            return reduction_statistic_message{
                .min_reduced_local_nodes =
                    std::min(s1.min_reduced_local_nodes, s2.min_reduced_local_nodes),
                .max_reduced_local_nodes =
                    std::max(s1.max_reduced_local_nodes, s2.max_reduced_local_nodes),
                .min_reduced_ghosts = std::min(s1.min_reduced_ghosts, s2.min_reduced_ghosts),
                .max_reduced_ghosts = std::max(s1.max_reduced_ghosts, s2.max_reduced_ghosts)};
          },
          ops::commutative));

  if (comm_.is_root()) {
    auto reduction_stats =
        std::vector<KernelStatistic::ReductionStatistic>(status_.reductions.size() + 1);
    // reduction-wise aggregate statistic
    for (std::size_t i = 0; i < status_.reductions.size(); ++i) {
      auto &stat = agg_reduction_statistics[i];
      auto &reduction = *status_.reductions[i];

      reduction_stats[i] = {.type = reduction.get_reduction_type(),
                            .min_reduced_local_nodes = stat.min_reduced_local_nodes,
                            .max_reduced_local_nodes = stat.max_reduced_local_nodes,
                            .min_reduced_ghosts = stat.min_reduced_ghosts,
                            .max_reduced_ghosts = stat.max_reduced_ghosts};
    }
    auto &received_stat = agg_reduction_statistics[status_.reductions.size()];
    reduction_stats[status_.reductions.size()] = {
        .type = ReductionType::received_dummy,
        .min_reduced_local_nodes = received_stat.min_reduced_local_nodes,
        .max_reduced_local_nodes = received_stat.max_reduced_local_nodes,
        .min_reduced_ghosts = received_stat.min_reduced_ghosts,
        .max_reduced_ghosts = received_stat.max_reduced_ghosts};
    auto &stat = agg_stat.value();
    return KernelStatistic{
        .kernel_graph = kernel_graph_stat.value(),
        .offset = stat.offset,
        .allreduce = stat.allreduce,
        .alltoallv = stat.alltoallv,
        .max_recv_border_updates = stat.max_recv_border_updates,
        .sum_recv_border_updates = stat.sum_recv_border_updates,
        .sum_recv_redundant_border_updates = stat.sum_recv_redundant_border_updates,
        .max_recv_redundant_border_updates_rel = stat.max_recv_redundant_border_updates_rel,
        .reduction_stats = std::move(reduction_stats)};
  } else {
    return std::nullopt;
  }
}

void Reducer::exclude_node(NodeID node) {
  // exclude
  KASSERT(status_.solution.node_status[node] == wis_status::UNSET ||
          status_.solution.node_status[node] == wis_status::DETACHED);
  status_.modified_queue.push_back(node);
  status_.solution.node_status[node] = wis_status::EXCLUDED;

  status_.graph.hide_node(node);
  if (status_.graph.is_active_border_node(node)) {
    on_border_status_change(node, wis_status::EXCLUDED);
  }
  if (!status_.graph.is_ghost(node)) {
    status_.reduction_measurements[active_reduction_index].reduced_local_nodes += 1;
  } else {
    status_.reduction_measurements[active_reduction_index].reduced_ghosts += 1;
  }
  // mark neighbors
  mark_neigh(node);
  if (status_.graph.is_ghost(node)) {
    for (auto border_vertex : status_.graph[node]) {
      if (!status_.graph.is_active_border_node(border_vertex)) {
        KASSERT(std::none_of(status_.graph[border_vertex].begin(),
                             status_.graph[border_vertex].end(),
                             [&](NodeID target) { return status_.graph.is_ghost(target); }));
        // border_vertex is not a border vertex anymore
        // this activates some reductions
        mark_neigh(border_vertex);
      }
    }
  }
}
void Reducer::fold_node(NodeID node) {
  // exclude
  KASSERT(status_.solution.node_status[node] == wis_status::UNSET ||
          status_.solution.node_status[node] == wis_status::DETACHED);
  status_.modified_queue.push_back(node);
  status_.solution.node_status[node] = wis_status::FOLDED;

  status_.graph.hide_node(node);
  if (!status_.graph.is_ghost(node)) {
    status_.reduction_measurements[active_reduction_index].reduced_local_nodes += 1;
  } else {
    status_.reduction_measurements[active_reduction_index].reduced_ghosts += 1;
  }

  // mark neighbors
  mark_neigh(node);
}
void Reducer::detach_node(NodeID node) {
  // exclude
  KASSERT(status_.solution.node_status[node] == wis_status::UNSET);
  KASSERT(status_.graph.is_active_border_node(node));
  status_.graph.hide_node(node);
  status_.reduction_measurements[active_reduction_index].reduced_local_nodes +=
      1;  // account at original rank

  status_.modified_queue.push_back(node);
  status_.solution.node_status[node] = wis_status::DETACHED;
  on_border_status_change(node, wis_status::DETACHED);
}
void Reducer::include_node(NodeID node) {
  // include node
  KASSERT(status_.solution.node_status[node] == wis_status::UNSET ||
          status_.solution.node_status[node] == wis_status::DETACHED);
  status_.modified_queue.push_back(node);
  status_.solution.node_status[node] = wis_status::INCLUDED;
  status_.graph.hide_node(node);
  if (status_.graph.is_active_border_node(node)) {
    on_border_status_change(node, wis_status::INCLUDED);
  }
  if (!status_.graph.is_ghost(node)) {
    status_.reduction_measurements[active_reduction_index].reduced_local_nodes += 1;
    status_.solution.solution_weight += status_.graph.get_weight(node);
  } else {
    status_.reduction_measurements[active_reduction_index].reduced_ghosts += 1;
  }
  // exclude neighbors
  for (auto &target : status_.graph[node]) {
    exclude_node(target);
  }
}
void Reducer::set_new_weight(NodeID node, NodeWeight new_weight) {
  status_.graph.set_weight(node, new_weight);
  if (status_.graph.is_active_border_node(node)) {
    on_border_weight_shift(node, new_weight);
  }
  mark_neigh(node);
}
void Reducer::set(NodeID node, wis_status new_node_status, bool poll) {
  if (new_node_status == wis_status::INCLUDED) {
    include_node(node);
  } else if (new_node_status == wis_status::FOLDED) {
    fold_node(node);
  } else if (new_node_status == wis_status::DETACHED) {
    detach_node(node);
  } else {
    exclude_node(node);
  }
  if (poll) {
    poll_updates_during_reduce_phase();
  }
}

void Reducer::include_node_with_bulk_hide(NodeID node, fast_set &neighbors, bool poll) {
  auto &graph = status_.graph;

  KASSERT(status_.solution.node_status[node] == wis_status::UNSET);
  status_.modified_queue.push_back(node);
  status_.solution.node_status[node] = wis_status::INCLUDED;
  graph.hide_node(node);
  if (graph.is_active_border_node(node)) {
    on_border_status_change(node, wis_status::INCLUDED);
  }
  if (!graph.is_ghost(node)) {
    status_.reduction_measurements[active_reduction_index].reduced_local_nodes += 1;
    status_.solution.solution_weight += graph.get_weight(node);
  } else {
    status_.reduction_measurements[active_reduction_index].reduced_ghosts += 1;
  }
  bulk_exclude(std::span(graph[node]), neighbors, poll);
}

void Reducer::bulk_exclude(std::span<const NodeID> exclude, fast_set &exclude_set, bool poll) {
  auto &graph = status_.graph;
  graph.bulk_hide(exclude, exclude_set);

  for (auto &node : exclude) {
    KASSERT(status_.solution.node_status[node] == wis_status::UNSET);
    status_.solution.node_status[node] = wis_status::EXCLUDED;
    status_.modified_queue.push_back(node);
    if (graph.is_active_border_node(node)) {
      on_border_status_change(node, wis_status::EXCLUDED);
    }
    if (!status_.graph.is_ghost(node)) {
      status_.reduction_measurements[active_reduction_index].reduced_local_nodes += 1;
    } else {
      status_.reduction_measurements[active_reduction_index].reduced_ghosts += 1;
    }
  }

  for (auto &modified : graph.get_last_bulk_modified_nodes()) {
    if (!graph.is_ghost(modified)) {
      mark_node(modified);
    }
  }

  if (poll) {
    poll_updates_during_reduce_phase();
  }
}

Reducer::ReduceStatus::ReduceStatus(Reducer::Config cfg, dist_dynamic_graph g,
                                    kamping::measurements::Timer<> *t)
    : block_rank(g.blocks()),
      solution(g.num_nodes(), g.num_ghosts()),
      graph(std::move(g)),
      cfg(cfg),
      reduce_timer(t),
      // modified_queue(graph.num_nodes() + graph.num_ghosts()),
      folded_queue(graph.num_nodes() + graph.num_ghosts()),
      alltoallv(0),
      allreduce(0),
      recv_border_updates(0),
      recv_redundant_border_updates(0) {
  std::iota(block_rank.begin(), block_rank.end(), 0);
  modified_queue.reserve(2 * (graph.num_nodes() + graph.num_ghosts()));
}

}  // namespace kadisredu::mwis
