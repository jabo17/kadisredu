//
// Created by jannickb on 11/27/24.
//

#pragma once

#include <kadisredu/algorithms/mwis/mwis_def.h>

namespace kadisredu::cli {

using namespace kadisredu;

/**
 * @brief parse CLI for KaDisRedu/mwis
 * @return -1 (success) if no error occurred during parsing the parameters; otherwise an error code
 * >= 0 is returned
 */
int parse_parameters(int argc, char *argv[], mwis::application_context &ctx, BlockID comm_size);

}  // namespace kadisredu::cli
