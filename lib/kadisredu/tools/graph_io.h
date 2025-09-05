//
// Created by jannickb on 1/6/25.
//

#pragma once

#include "kadisredu/definitions.h"

#include <mpi.h>
#include <kagen/kagen.h>

namespace kadisredu {
void write_metis_graph(MPI_Comm mpi_comm, kagen::Graph &graph, std::string filename);
}