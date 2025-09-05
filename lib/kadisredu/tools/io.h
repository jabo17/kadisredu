#pragma once

#include <fstream>
#include <kamping/communicator.hpp>

namespace kadisredu {
// Based on dist_io.h from dKaMinPar by D. Seemaier
template <typename T>
void write_distributed_vector(const kamping::Communicator<> &comm, const std::string &filename,
                              const std::vector<T> &vec) {
  using namespace kamping;

  if (comm.is_root()) {
    std::ofstream out(filename, std::ios::trunc);
  }

  for (int pe = 0; pe < comm.size(); ++pe) {
    if (pe == comm.rank()) {
      std::ofstream out(filename, std::ios::app);
      for (const auto &val : vec) {
        out << val << std::endl;
      }
    }
    comm.barrier();
  }
}
}  // namespace kadisredu
