
// Created by Jannick Borowitz on 29.05.24.
//

#include "cli.h"

#include <fstream>
#include <kadisredu/algorithms/mwis/kadisredu_mwis.h>
#include <kadisredu/algorithms/mwis/mwis_def.h>
#include <kadisredu/graphutils/dist_dynamic_graph_factories.h>
#include <kadisredu/tools/io.h>
#include <kadisredu/tools/json_printer.h>
#include <kadisredu/tools/logger.h>
#include <kagen/io.h>
#include <kagen/tools/converter.h>
#include <kamping/checking_casts.hpp>
#include <kamping/communicator.hpp>
#include <kamping/environment.hpp>
#include <kamping/measurements/printer.hpp>
#include <mpi.h>
#include <nlohmann/json.hpp>

namespace {
kagen::Graph load_kagen_graph(const kadisredu::mwis::application_context &app) {
  using namespace kagen;

  KaGen generator(MPI_COMM_WORLD);
  generator.UseCSRRepresentation();
  if (app.check_input_graph) {
    generator.EnableUndirectedGraphVerification();
  }
  if (app.experiment) {
    generator.EnableBasicStatistics();
    generator.EnableOutput(true);
  }

  return [&] {
    if (std::find(app.kagen_option_string.begin(), app.kagen_option_string.end(), ';') !=
        app.kagen_option_string.end()) {
      return generator.GenerateFromOptionString(app.kagen_option_string);
    } else {
      return generator.ReadFromFile(app.kagen_option_string, FileFormat::EXTENSION,
                                    GraphDistribution::BALANCE_VERTICES);
    }
  }();
}

bool check_solution(MPI_Comm mpi_comm, kadisredu::GlobalNodeWeight local_wis_sol,
                    const std::vector<kamping::kabool> &node_status, const kagen::Graph &graph) {
  using namespace kagen;
  using namespace kaminpar;
  using namespace kaminpar::dist;
  using namespace kamping;
  KASSERT(node_status.size() == graph.NumberOfLocalVertices());

  KADISREDU_RLOG << "Verifying solution ...";

  // do something with kamping
  kamping::Communicator<> comm(mpi_comm);

  auto ghost_mapping = kadisredu::build_ghost_mapping(comm, graph);

  // get status for ghosts of KaGen graph
  struct ghost_data {
    GlobalNodeID ghost;
    bool status;
  };

  std::vector<std::vector<ghost_data>> send_border_data(comm.size());
  // prepare one message for each border node and adjacent block
  kadisredu::fast_set reached_blocks(comm.size());
  auto n = graph.NumberOfLocalVertices();
  auto [begin_v, end_v] = graph.vertex_range;

  for (NodeID v = 0; v < n; v++) {
    reached_blocks.clear();
    for (kagen::SInt e = graph.xadj[v]; e < graph.xadj[v + 1]; ++e) {
      auto u = static_cast<GlobalNodeID>(graph.adjncy[e]);
      if (u < begin_v || u >= end_v) {
        // u is a ghost
        auto reached_ghost = ghost_mapping.get(u);
        if (reached_blocks.add(ghost_mapping.blocks[reached_ghost])) {
          send_border_data[ghost_mapping.blocks[reached_ghost]].push_back(
              {.ghost = v + begin_v, .status = node_status[v]});
        }
      }
    }
  }

  // we expect from every ghost a message
  std::vector<int> recv_ghost_message_counts(comm.size(), 0);
  for (auto &rank : ghost_mapping.blocks) {
    recv_ghost_message_counts[rank]++;
  }
  auto recv_ghost_data = with_flattened(send_border_data, comm.size()).call([&](auto... flattened) {
    return comm.alltoallv(std::move(flattened)..., recv_counts(recv_ghost_message_counts));
  });
  KASSERT(recv_ghost_data.size() == ghost_mapping.num_ghosts());

  // set ghost status
  std::vector<bool> ghosts_status(ghost_mapping.num_ghosts());
  for (auto [global_node, ghost_status] : recv_ghost_data) {
    ghosts_status[ghost_mapping.get(global_node)] = ghost_status;
  }

  // check maximal weight independent set locally
  GlobalNodeWeight local_check_wis_weight = 0;
  bool check_failed = false;
  for (NodeID v = 0; v < n; v++) {
    if (node_status[v]) {
      local_check_wis_weight += static_cast<GlobalNodeWeight>(graph.vertex_weights[v]);
      // v is included, check that no neighbor exists in the solution
      bool independent = true;
      for (kagen::SInt edge = graph.xadj[v]; edge < graph.xadj[v + 1]; edge++) {
        auto u = static_cast<GlobalNodeID>(graph.adjncy[edge]);
        if (u < graph.vertex_range.first || u >= graph.vertex_range.second) {
          if (ghosts_status[ghost_mapping.get(u)]) {
            independent = false;
            break;
          }
        } else {
          auto local_u = u - graph.vertex_range.first;
          if (node_status[local_u]) {
            independent = false;
            break;
          }
        }
      }
      if (!independent) {
        KADISREDU_LOG_ERROR << "Independent set check failed!";
        check_failed = true;
        break;
      }
    } else {
      // v is not in the solution; ensure it is no free node
      bool free_node = true;
      bool border_node = false;
      for (kagen::SInt edge = graph.xadj[v]; edge < graph.xadj[v + 1]; edge++) {
        auto u = static_cast<GlobalNodeID>(graph.adjncy[edge]);
        if (u < graph.vertex_range.first || u >= graph.vertex_range.second) {
          border_node = true;
          if (ghosts_status[ghost_mapping.get(u)]) {
            free_node = false;
            break;
          }
        } else {
          auto local_u = u - graph.vertex_range.first;
          if (node_status[local_u]) {
            free_node = false;
            break;
          }
        }
      }
      if (free_node) {
        KADISREDU_LOG_ERROR << "Maximality check failed! (border vertex?: " << border_node << ") "
                            << (graph.xadj[v + 1] - graph.xadj[v]);
        check_failed = true;
        break;
      }
    }
  }

  if (comm.allreduce_single(send_buf((kabool)check_failed), op(ops::logical_or<>()))) {
    return false;
  }

  auto compare_weights =
      comm.allreduce(send_buf({local_check_wis_weight, local_wis_sol}), op(ops::plus<>()));

  if (compare_weights[0] != compare_weights[1]) {
    if (comm.is_root()) {
      KADISREDU_LOG_ERROR << "The overall solution weight does not match the summed weights of the "
                             "provided solution: "
                          << compare_weights[0] << " (sum) != " << compare_weights[1]
                          << " (provided sum)";
    }
    return false;
  } else {
    KADISREDU_RLOG << "Verified WIS weight: " << compare_weights[0];
  }

  return true;
}

void optimize_reduction_phases(kadisredu::mwis::application_context &app) {
  // Optimize reduction phase order by removing duplicate reduction phases
  // that do not partition in-between.
  auto new_red_ctx_size = app.red_ctx.size();
  std::size_t i = 1;
  while (i < new_red_ctx_size) {
    auto &first = app.red_ctx[i - 1];
    auto &second = app.red_ctx[i];

    if (first == second && second.initial_partitioning.disable_partitioning) {
      using std::swap;
      --new_red_ctx_size;
      swap(second, app.red_ctx[new_red_ctx_size]);
    } else {
      ++i;
    }
  }
  if (new_red_ctx_size != app.red_ctx.size()) {
    KADISREDU_RLOG << "Removed " << app.red_ctx.size() - new_red_ctx_size
                   << " redundant reduction phase(s).";
    app.red_ctx.resize(new_red_ctx_size);
  }
}

void disable_partitioning(kadisredu::mwis::application_context &app) {
  for (auto &red_ctx : app.red_ctx) {
    red_ctx.initial_partitioning.disable_partitioning = true;
  }
  app.greedy_ctx.initial_partitioning.disable_partitioning = true;
}

void write_json_stats_on_root(const kamping::Communicator<> &comm,
                              kadisredu::mwis::application_context &app,
                              kadisredu::mwis::kadisredu_mwis &solver,
                              kamping::measurements::Timer<> &t) {
  using namespace kadisredu;
  using namespace kamping;

  // timer-tree stats
  std::stringstream timer_sjson;
  t.aggregate_and_print(measurements::SimpleJsonPrinter{timer_sjson});

  if (comm.is_root()) {
    std::ofstream out;
    auto &filename = app.json_output_path;
    if (!filename.empty() && comm.is_root()) {
      out.open(filename);
      if (!out.is_open()) {
        std::cerr << "Failed to open json_output_path " << filename << std::endl;
        comm.abort();
      }
    }

    nlohmann::json json_dump;
    json_dump["solver"] = solver.take_json_dump();
    json_dump["timer"] = JsonPrinter::get_json_from_stream(timer_sjson)["data"]["root"];
    JsonPrinter::pretty_print_json(json_dump, out);

    // close file stream
    if (out.is_open()) {
      out.close();
    }
  }
}

void warmup_mpi(MPI_Comm mpi_comm) {
  using namespace kamping;
  kamping::Communicator comm(mpi_comm);

  KADISREDU_RLOG << "Running MPI-warmup";
  std::vector<int> send;
  std::vector<int> counts(comm.size(), 0);

  for (int k = 0; k < 2; k++) {
    for (int i = 0; i < comm.size(); i++) {
      for (int j = 0; j < (k + i + comm.rank()) % comm.size(); j++) {
        counts[i]++;
        send.push_back(j);
      }
    }
    auto recv = comm.alltoallv(send_buf(send), send_counts(counts));
    asm volatile("" : : "r,m"(recv) : "memory");
  }
  KADISREDU_RLOG << "Finished MPI-warmup";
}

}  // namespace

int main(int argc, char *argv[]) {
  using namespace kagen;
  using namespace kaminpar;
  using namespace kaminpar::dist;
  using namespace kamping;
  using namespace kadisredu;
  using namespace kadisredu::mwis;

  // do something with kamping
  kamping::Environment e;
  kamping::Communicator comm;

  // setup CLI
  application_context app;

  int err = cli::parse_parameters(argc, argv, app, comm.size());
  if (err != -1) {
    return 1;
  }

  if (comm.size() == 1) {
    disable_partitioning(app);
    KADISREDU_RLOG << "Turned off partitioner since #cores=1.";
  }
  optimize_reduction_phases(app);

  if (app.warmup_mpi) {
    warmup_mpi(comm.mpi_communicator());
  }

  timer t_elapsed;
  t_elapsed.restart();
  auto &t = measurements::timer();  // singleton instance

  // ##################################################################
  KADISREDU_RLOG << "Loading graph ...";
  t.start("KaGen");
  // Load the graph via KaGen
  auto graph = load_kagen_graph(app);
  t.stop();
  KADISREDU_RLOG << "Loaded graph";

  // ##################################################################
  KADISREDU_RLOG << "Running KaDisRedu ...";
  t.start("KaDisRedu");

  // run solver scheme
  auto solver = kadisredu_mwis(comm.mpi_communicator(), app);
  auto n = graph.NumberOfLocalVertices();
  auto solution = solver.reduce_and_greedy(std::move(graph));
  KASSERT(solution.node_status.size() == n);
  // ##################################################################
  KADISREDU_RLOG << "Setting solution ...";
  t.start("set-solution");

  // check solution for KaGen graph
  // obtain solution status for original KaGen graph
  std::vector<kabool> independent_set(solution.node_status.size());
  std::ranges::transform(solution.node_status, independent_set.begin(),
                         [](wis_status status) { return status == wis_status::INCLUDED; });
  t.stop();  // stop set-solution
  t.stop();  // stop KaDisRedu
  KADISREDU_RLOG << "Finished KaDisRedu!";
  // ##################################################################
  if (app.print_independent_set) {
    write_distributed_vector(comm, app.independent_set_filename, independent_set);
  }
  // ##################################################################
  KADISREDU_RLOG << "Writing statistics ...";
  write_json_stats_on_root(comm, app, solver, t);
  // ##################################################################
  KADISREDU_RLOG << "Checking solution";
  // checking solution on input kagen graph
  graph = load_kagen_graph(app);  // re-load graph (might have been freed)
  if (!check_solution(comm.mpi_communicator(), solution.solution_weight, independent_set, graph) &&
      comm.root()) {
    KADISREDU_LOG_ERROR << "Solution is wrong (see log above)";
    return 0;
  }
  KADISREDU_RLOG << "Verified solution! ";

  return 0;
}
