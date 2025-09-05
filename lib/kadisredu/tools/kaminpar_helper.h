//
// Created by borowitz on 30.05.25.
//

#pragma once

#include "kadisredu/algorithms/mwis/mwis_def.h"

#include <dkaminpar.h>
#include <kagen/kagen.h>
#include <kamping/communicator.hpp>

namespace kadisredu {
class kaminpar_helper {
 private:
  static void import_kagen_graph_into_partitioner(kagen::Graph graph,
                                                  kaminpar::dKaMinPar &partitioner) {
    using namespace kagen;
    using namespace kaminpar;
    using namespace kaminpar::dist;

    // We use `unsigned long` here since we currently do not have any MPI type
    // definitions for GlobalNodeID
    static_assert(std::is_same_v<GlobalNodeID, unsigned long>);
    std::vector<GlobalNodeID> vtxdist =
        BuildVertexDistribution<unsigned long>(graph, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    // ... if the data types are not the same, we would need to re-allocate memory
    // for the graph; to this if we ever need it ...
    std::vector<SInt> xadj = graph.TakeXadj<>();
    std::vector<SInt> adjncy = graph.TakeAdjncy<>();
    // std::vector<SSInt> vwgt = graph.TakeVertexWeights<>();
    std::vector<SSInt> adjwgt = graph.TakeEdgeWeights<>();

    static_assert(sizeof(SInt) == sizeof(GlobalNodeID));
    static_assert(sizeof(SInt) == sizeof(GlobalEdgeID));
    static_assert(sizeof(SSInt) == sizeof(GlobalNodeWeight));
    static_assert(sizeof(SSInt) == sizeof(GlobalEdgeWeight));

    partitioner.copy_graph(vtxdist, {reinterpret_cast<GlobalEdgeID *>(xadj.data()), xadj.size()},
                           {reinterpret_cast<GlobalNodeID *>(adjncy.data()), adjncy.size()}, {},
                           {reinterpret_cast<GlobalEdgeWeight *>(adjwgt.data()), adjwgt.size()});
  }

 public:
  static void partition_graph(kagen::Graph graph, MPI_Comm comm, int seed,
                              const kadisredu::mwis::partitioner_context &part_ctx,
                              std::vector<kaminpar::dist::BlockID> &partition) {
    using namespace kaminpar;
    using namespace kaminpar::dist;
    using namespace kamping;

    Communicator<> k_comm(comm);
    GlobalNodeID g_n = graph.NumberOfGlobalVertices();
    if (g_n < part_ctx.minimum_vertices_per_rank * k_comm.size()) {
      auto n = graph.vertex_range.second - graph.vertex_range.first;
      partition.resize(n);
      for (NodeID v = 0; v < n; ++v) {
        partition[v] = k_comm.root();
      }
      return;
    }

    // partitioning context
    Context ctx = create_default_context();
    // Apply partitioner in non-trivial case (k>1)
    dKaMinPar partitioner(comm, 1, ctx);
    partitioner.set_output_level(OutputLevel::QUIET);
    dKaMinPar::reseed(seed);
    // Import distributed graph
    import_kagen_graph_into_partitioner(std::move(graph), partitioner);
    partitioner.compute_partition(static_cast<PEID>(k_comm.size()), part_ctx.epsilon, partition);
  }
};
}  // namespace kadisredu
