//
// Created by jannickb on 1/6/25.
//

#include "kadisredu/tools/graph_io.h"
#include "kadisredu/tools/logger.h"

#include <kamping/communicator.hpp>
#include <kagen/tools/converter.h>
#include <kagen/io/metis.h>

#include <fstream>

namespace kadisredu {
void write_metis_graph(MPI_Comm mpi_comm, kagen::Graph &graph, std::string filename) {
  using namespace kamping;
  using namespace kagen;

  kamping::Communicator<> comm(mpi_comm);

  // print static graph to file
  using namespace kagen;

  std::vector<GlobalNodeID> vtxdist =
      BuildVertexDistribution<unsigned long>(graph, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  OutputGraphConfig out_cfg{};
  out_cfg.filename = filename;
  out_cfg.formats[0] = FileFormat::METIS;
  GraphInfo info{graph, comm.mpi_communicator()};

  graph.edges = BuildEdgeListFromCSR(graph.vertex_range, graph.xadj, graph.adjncy);
  MetisWriter writer{out_cfg, graph, info, static_cast<PEID>(comm.rank()),
                     static_cast<PEID>(comm.size())};
  WriteGraph(writer, out_cfg, true, comm.mpi_communicator());
  KADISREDU_RLOG << "Written graph to " << filename;
}
}