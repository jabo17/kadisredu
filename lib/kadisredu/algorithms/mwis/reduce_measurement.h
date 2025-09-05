//
// Created by jannickb on 10/19/24.
//

#pragma once

#include "kadisredu/definitions.h"
#include "mwis_def.h"

namespace kadisredu::mwis {

class ReduceMeasurement {
 public:
  NodeID reduced_local_nodes = 0;
  NodeID reduced_ghosts = 0;
};

}  // namespace kadisredu::mwis