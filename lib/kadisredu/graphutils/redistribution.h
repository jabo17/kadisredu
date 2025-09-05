#pragma once

#include "kadisredu/data_structures/dist_dynamic_graph.h"
#include "kadisredu/definitions.h"
#include "kadisredu/tools/logger.h"

#include <algorithm>
#include <kamping/communicator.hpp>

namespace kadisredu::redistribution {

struct ghost_node_mapper {
  ghost_node_mapper() = default;
  explicit ghost_node_mapper(std::vector<GlobalNodeID> node_dist, NodeID reserve)
      : map(reserve), node_dist(std::move(node_dist)) {}

  dist_dynamic_graph::ghost_node_mapping map;
  std::vector<GlobalNodeID> node_dist;

  NodeID add(GlobalNodeID node) {
    auto [ghost, added] = map.add(node, 0);
    if (added) {
      map.blocks[ghost] = std::distance(node_dist.begin(),
                                        std::upper_bound(node_dist.begin(), node_dist.end(), node));
      KASSERT(node < node_dist[map.blocks[ghost]]);
      KASSERT(map.blocks[ghost] == 0 || node >= node_dist[map.blocks[ghost] - 1]);
      KASSERT(kamping::comm_world().size() != 2 ||
              !static_cast<bool>(kamping::comm_world().rank()) == map.blocks[ghost]);
    }
    return ghost;
  }
};

struct RedistributionRestoreData {
  std::vector<NodeID> new_to_old;
  std::vector<int> partition_counts;
  std::vector<int> reverse_partition_counts;
  std::vector<NodeWeight> old_weights;

  RedistributionRestoreData() = default;
  RedistributionRestoreData(NodeID local_vertices, BlockID partitions)
      : new_to_old(local_vertices), partition_counts(partitions, 0) {}
};

struct VertexRange {
  GlobalNodeID begin;
  GlobalNodeID end;
};

struct EdgeData {
  NodeID source;
  NodeID target;
  BlockID target_block;
};

struct GraphData {
  std::vector<NodeWeight> weights;
  std::vector<EdgeData> edges;
  VertexRange vertex_range;
  std::vector<GlobalNodeID> global_partition_displ;

  void free_edges() { free_vec(edges); }

  void free_weights() { free_vec(weights); }

 private:
  template <typename T>
  void free_vec(std::vector<T> &vec) {
    std::vector<T>().swap(vec);
  }
};

inline dist_dynamic_graph build_dist_dynamic_graph(const kamping::Communicator<> &comm,
                                                   GraphData data) {
  KADISREDU_RLOG << "Build redistributed dynamic graph";

  auto [begin_v, end_v] = data.vertex_range;
  NodeID n = end_v - begin_v;

  auto global_id = [&](NodeID rank_in_part, BlockID block) {
    return data.global_partition_displ[block] + rank_in_part;
  };

  // build ghost mapping
  dist_dynamic_graph::ghost_node_mapping ghosts;
  for (auto &edge : data.edges) {
    if (edge.target_block != comm.rank()) {
      GlobalNodeID g_target = global_id(edge.target, edge.target_block);
      auto [ghost, added] = ghosts.add(g_target, edge.target_block);
      edge.target = ghost + n;
    } else {
      KASSERT(edge.target <= n);
    }
  }
  KADISREDU_RLOG << "[build dyn graph] " << "built ghost mapping";

  // init graph
  dist_dynamic_graph graph(comm.mpi_communicator(), begin_v, n, std::move(ghosts));

  // insert edges
  auto &edges = data.edges;
  auto read_pos = edges.begin();
  while (read_pos != edges.end()) {
    NodeID source = read_pos->source;
    auto end_source_edges = std::find_if(read_pos, edges.end(),
                                         [&source](auto &edge) { return edge.source != source; });

    std::vector<NodeID> neighbors_buffer(std::distance(read_pos, end_source_edges));
    std::transform(read_pos, end_source_edges, neighbors_buffer.begin(),
                   [](const auto &edge) { return edge.target; });
    graph.build_neighborhood(source, std::move(neighbors_buffer));
    read_pos = end_source_edges;
  }
  KADISREDU_RLOG << "[build dyn graph] " << "built neighborhoods";
  data.free_edges();

  // set weights
  for (NodeID v = 0; v < n; v++) {
    graph.set_weight(v, data.weights[v]);
  }
  data.free_weights();

  // build border
  graph.build_border();
  KADISREDU_RLOG << "[build dyn graph] " << "built ghosts' neighborhoods";

  // sync ghosts
  auto ghost_weights = graph.sync_ghosts_with_lazy_updates<NodeWeight>(
      graph.get_border(),
      [&](NodeID v) -> std::optional<NodeWeight> { return graph.get_weight(v); });
  std::ranges::for_each(
      ghost_weights, [&](auto &up) { graph.set_weight(graph.get_ghost(up.global_id), up.update); });
  KADISREDU_RLOG << "[build dyn graph] " << "finished setting vertex weights of ghosts";

  KADISREDU_RLOG << "[build dyn graph] " << "Finished";
  return graph;
}

inline RedistributionRestoreData get_redistribution_restore_data(
    const kamping::Communicator<> &comm, VertexRange vertex_range,
    const std::vector<BlockID> &partition) {
  using namespace kamping;
  auto [begin_v, end_v] = vertex_range;
  NodeID n = end_v - begin_v;

  RedistributionRestoreData restore_data(n, comm.size());

  auto &partition_counts = restore_data.partition_counts;
  std::for_each(partition.begin(), partition.end(), [&](BlockID b) { partition_counts[b]++; });

  std::vector<NodeID> partition_displs(comm.size(), 0);
  std::exclusive_scan(partition_counts.begin(), partition_counts.end(), partition_displs.begin(),
                      0);

  auto &new_to_old = restore_data.new_to_old;  // local vertices are ordered by partition assignment
  {
    auto partition_displs_copy = partition_displs;
    for (NodeID v = 0; v < n; ++v) {
      NodeID new_v = partition_displs_copy[partition[v]]++;
      new_to_old[new_v] = v;
    }
  }

  // redistribute vertices
  auto &reverse_partition_counts = restore_data.reverse_partition_counts;
  reverse_partition_counts = comm.alltoall(send_buf(partition_counts));

  return restore_data;
}

inline auto redistribute_graph(const kamping::Communicator<> &comm, kagen::Graph &&graph,
                               const std::vector<BlockID> &partition) {
  KASSERT(graph.NumberOfLocalVertices() == partition.size());

  using namespace kamping;

  NodeID n = graph.NumberOfLocalVertices();
  auto [begin_v, end_v] = graph.vertex_range;

  RedistributionRestoreData restore_data =
      get_redistribution_restore_data(comm, {.begin = begin_v, .end = end_v}, partition);
  auto &partition_counts = restore_data.partition_counts;
  auto &reverse_partition_counts = restore_data.reverse_partition_counts;
  auto &new_to_old = restore_data.new_to_old;

  GraphData graph_data;

  std::vector<NodeID> old_to_new(n);  // reverse mapping
  {
    for (NodeID v = 0; v < n; ++v) {
      old_to_new[new_to_old[v]] = v;
    }
  }

  // compute rank of vertices in partition
  NodeID recv_n = std::reduce(reverse_partition_counts.begin(), reverse_partition_counts.end(), 0);
  std::vector<NodeID> id(recv_n);
  std::iota(id.begin(), id.end(), 0);
  auto rank_in_partition = comm.alltoallv(send_buf(id), send_counts(reverse_partition_counts),
                                          recv_counts(partition_counts));

  // compute new vertex distribution
  std::vector<GlobalNodeID> global_partition_counts(comm.size(), 0);
  KASSERT(partition_counts.size() == comm.size());
  KASSERT(std::reduce(partition_counts.begin(), partition_counts.end(), 0) == n);
  std::ranges::transform(partition_counts, global_partition_counts.begin(),
                         [](int x) -> GlobalNodeID { return static_cast<GlobalNodeID>(x); });
  comm.allreduce_inplace(send_recv_buf(global_partition_counts), op(ops::plus<>()));
  KASSERT(std::reduce(global_partition_counts.begin(), global_partition_counts.end(),
                      GlobalNodeID{0}) == graph.NumberOfGlobalVertices());

  auto &global_partition_displ = graph_data.global_partition_displ;
  global_partition_displ = global_partition_counts;
  std::exclusive_scan(global_partition_counts.begin(), global_partition_counts.end(),
                      global_partition_displ.begin(), 0);
  // recv global vertex range
  GlobalNodeID recv_begin_vtx = global_partition_displ[comm.rank()];
  GlobalNodeID recv_end_vtx = recv_begin_vtx + global_partition_counts[comm.rank()];
  graph_data.vertex_range = {.begin = recv_begin_vtx, .end = recv_end_vtx};

  // redistribute weights
  std::vector<NodeWeight> vertex_weights(n);
  for (NodeID v = 0; v < n; ++v) {
    // rearrange weights locally
    vertex_weights[old_to_new[v]] = static_cast<NodeWeight>(graph.vertex_weights[v]);
  }
  graph_data.weights = comm.alltoallv(send_buf(vertex_weights), send_counts(partition_counts),
                                      recv_counts(reverse_partition_counts));

  // partition of ghosts
  const GlobalNodeID graph_ghost_count_estimation = n;
  // vertex distribution for ghosts mapper stores displacement + count for each partition
  std::vector<GlobalNodeID> vertex_distribution(comm.size());
  vertex_distribution[comm.rank()] = end_v;
  comm.allgather_inplace(send_recv_buf(vertex_distribution));
  ghost_node_mapper graph_ghosts(std::move(vertex_distribution), graph_ghost_count_estimation);

  fast_set graph_interface_vertices(n);
  for (NodeID v = 0; v < n; ++v) {
    for (kagen::SInt e = graph.xadj[v]; e < graph.xadj[v + 1]; ++e) {
      auto u = static_cast<GlobalNodeID>(graph.adjncy[e]);
      if (u < begin_v || u >= end_v) {
        // ghost
        graph_ghosts.add(u);
        graph_interface_vertices.add(v);
      }
    }
  }

  struct PartitionMessage {
    GlobalNodeID id;
    NodeID rank_in_part;
    BlockID part;
  };
  std::vector<std::vector<PartitionMessage>> send_partition_buckets(comm.size());
  fast_set reached_block(comm.size());
  for (NodeID v = 0; v < n; ++v) {
    if (graph_interface_vertices.contains(v)) {
      reached_block.clear();
      BlockID reached_blocks = 0;
      for (kagen::SInt e = graph.xadj[v]; e < graph.xadj[v + 1]; ++e) {
        auto u = static_cast<GlobalNodeID>(graph.adjncy[e]);
        if (u < begin_v || u >= end_v) {
          auto b = graph_ghosts.map.get_block_from_global(u);
          if (reached_block.add(b)) {
            send_partition_buckets[b].push_back(
                {v + graph.vertex_range.first, rank_in_partition[old_to_new[v]], partition[v]});
          }
        }
      }
    }
  }
  std::vector<BlockID> ghost_partition(graph_ghosts.map.num_ghosts());
  std::vector<NodeID> ghost_rank_in_partition(graph_ghosts.map.num_ghosts());
  {
    auto graph_ghost_data = with_flattened(send_partition_buckets).call([&](auto... flattened) {
      return comm.alltoallv(std::move(flattened)...);
    });
    for (auto [id, rank_in_part, part] : graph_ghost_data) {
      auto lid = graph_ghosts.map.get(id);
      ghost_rank_in_partition[lid] = rank_in_part;
      ghost_partition[lid] = part;
    }
  }

  // edge counts per partition
  std::vector<int> edge_counts(comm.size(), 0);
  for (NodeID v = 0; v < n; v++) {
    edge_counts[partition[v]] += static_cast<int>(graph.xadj[v + 1] - graph.xadj[v]);
  }
  std::vector<EdgeID> edge_displ(comm.size(), 0);
  std::exclusive_scan(edge_counts.begin(), edge_counts.end(), edge_displ.begin(), 0);

  // edges
  std::vector<EdgeData> edges(graph.NumberOfLocalEdges());
  {
    auto copy_edge_displ = edge_displ;
    for (NodeID v = 0; v < n; v++) {
      BlockID p_source = partition[v];
      NodeID source = rank_in_partition[old_to_new[v]];
      for (kagen::SInt e = graph.xadj[v]; e < graph.xadj[v + 1]; e++) {
        auto u = static_cast<GlobalNodeID>(graph.adjncy[e]);

        NodeID target;
        BlockID p_target;
        if (u < begin_v || u >= end_v) {
          // u is a ghost
          auto lu = graph_ghosts.map.get(u);
          p_target = ghost_partition[lu];
          target = ghost_rank_in_partition[lu];
        } else {
          // u is local
          auto lu = u - begin_v;
          p_target = partition[lu];
          target = rank_in_partition[old_to_new[lu]];
        }
        KASSERT(target < global_partition_counts[p_target]);
        edges[copy_edge_displ[p_source]++] = {
            .source = source, .target = target, .target_block = p_target};
      }
    }
  }
  // consume graph
  restore_data.old_weights.resize(n);
  std::ranges::transform(graph.vertex_weights, restore_data.old_weights.begin(),
                         [](auto weight) { return static_cast<NodeWeight>(weight); });
  graph.Clear();
  graph_data.edges = comm.alltoallv(send_buf(edges), send_counts(edge_counts));

  return std::make_pair(restore_data, graph_data);
}

inline std::pair<RedistributionRestoreData, dist_dynamic_graph>
redistribute_and_build_dist_dynamic_graph(MPI_Comm communicator, kagen::Graph &&graph,
                                          const std::vector<BlockID> &partition) {
  kamping::Communicator<> comm(communicator);
  // Note that redistribute_graph consumes `graph'
  auto [restore_data, graph_data] = redistribute_graph(comm, std::move(graph), partition);
  return std::make_pair(restore_data, build_dist_dynamic_graph(comm, std::move(graph_data)));
}

}  // namespace kadisredu::redistribution
