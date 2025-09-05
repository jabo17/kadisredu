//
// Created by jannickb on 6/17/24.
//

#pragma once

#include "kadisredu/data_structures/dist_dynamic_graph.h"
#include "kadisredu/definitions.h"
#include "kadisredu/tools/logger.h"

#include <kamping/communicator.hpp>
#include <kamping/measurements/timer.hpp>
#include <mpi.h>
#include <numeric>

namespace kadisredu {

inline dist_dynamic_graph::ghost_node_mapping build_ghost_mapping(
    const kamping::Communicator<>& comm, const kagen::Graph& graph) {
  using namespace kamping;
  struct ghost_node_mapper {
    ghost_node_mapper() = default;
    explicit ghost_node_mapper(std::vector<GlobalNodeID> node_dist, NodeID reserve)
        : map(reserve), node_dist(std::move(node_dist)) {}

    dist_dynamic_graph::ghost_node_mapping map;
    std::vector<GlobalNodeID> node_dist;

    NodeID add(GlobalNodeID node) {
      auto [ghost, added] = map.add(node, 0);
      if (added) {
        map.blocks[ghost] = std::distance(
            node_dist.begin(), std::upper_bound(node_dist.begin(), node_dist.end(), node));
        KASSERT(node < node_dist[map.blocks[ghost]]);
        KASSERT(kamping::comm_world().size() != 2 ||
                !static_cast<bool>(kamping::comm_world().rank()) == map.blocks[ghost]);
      }
      return ghost;
    }
  };

  NodeID n = graph.NumberOfLocalVertices();
  auto [begin_v, end_v] = graph.vertex_range;

  // built vertex distribution for ghost mapper
  std::vector<GlobalNodeID> vertex_distribution(comm.size(), 0);
  vertex_distribution[comm.rank()] = end_v;
  comm.allgather_inplace(send_recv_buf(vertex_distribution));

  // init ghost mapper
  ghost_node_mapper local_ghosts(std::move(vertex_distribution), comm.size() == 1 ? 0 : n);

  if (comm.size() > 1) {
    for (NodeID v = 0; v < n; ++v) {
      for (kagen::SInt e = graph.xadj[v]; e < graph.xadj[v + 1]; ++e) {
        auto u = static_cast<GlobalNodeID>(graph.adjncy[e]);
        // check if u is a ghost
        if (u < begin_v || u >= end_v) {
          local_ghosts.add(u);
        }
      }
    }
  }

  return std::move(local_ghosts.map);
}

inline dist_dynamic_graph build_dist_dynamic_graph(const kamping::Communicator<>& comm,
                                                   const kagen::Graph& graph,
                                                   kamping::measurements::Timer<>& t) {
  using namespace kamping;

  auto n = graph.NumberOfLocalVertices();
  auto [begin_v, end_v] = graph.vertex_range;
  // insert ghosts
  t.start("ghost-mapper-insert");
  dist_dynamic_graph::ghost_node_mapping ghost_mapping = build_ghost_mapping(comm, graph);
  KADISREDU_RLOG << "[build dyn graph] " << "built vertex dist.";
  t.stop();
  KADISREDU_RLOG << "[build dyn graph] " << "built ghost mapping";

  dist_dynamic_graph extracted(comm.mpi_communicator(), begin_v, n, std::move(ghost_mapping));

  t.start("build-neighborhoods");
  // insert edges
  for (NodeID v = 0; v < n; v++) {
    const auto offset = graph.xadj[v];
    std::vector<NodeID> neighbors_buffer(graph.xadj[v + 1] - offset);
    for (kagen::SInt e = offset; e < graph.xadj[v + 1]; e++) {
      KASSERT(e < graph.adjncy.size());
      const auto u = static_cast<GlobalNodeID>(graph.adjncy[e]);
      // is u a ghost
      if (u < graph.vertex_range.first || u >= graph.vertex_range.second) {
        neighbors_buffer[e - offset] = extracted.get_ghost(u);  //
      } else {
        neighbors_buffer[e - offset] = u - begin_v;
      }
    }
    extracted.build_neighborhood(v, std::move(neighbors_buffer));
  }
  t.stop();

  KADISREDU_RLOG << "[build dyn graph] " << "built neighborhoods";

  for (NodeID v = 0; v < n; ++v) {
    extracted.set_weight(v, static_cast<NodeWeight>(graph.vertex_weights[v]));
  }

  t.start("build-border");
  extracted.build_border();
  t.stop();

  // sync ghosts
  t.start("exchange-weights-for-ghosts");
  auto ghost_updates = extracted.sync_ghosts_with_lazy_updates<NodeWeight>(
      extracted.get_border(), [&](NodeID v) { return extracted.get_weight(v); });
  t.stop();  // exchanged ghost weights
  t.start("set-ghost-weights");
  std::for_each(ghost_updates.begin(), ghost_updates.end(), [&extracted](auto& up) {
    extracted.set_weight(extracted.get_ghost(up.global_id), up.update);
  });
  t.stop();

  KADISREDU_RLOG << "[build dyn graph] " << "built ghosts' neighborhoods and set vertex weights";
  KADISREDU_RLOG << "Finished building dynamic graph!";

  return extracted;
}

}  // namespace kadisredu
