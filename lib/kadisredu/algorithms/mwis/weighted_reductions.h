//
// Created by jannickb on 6/24/24.
//

#pragma once

#include "kadisredu/algorithms/mwis/mwis_def.h"
#include "kadisredu/algorithms/mwis/reduction_ptr.h"
#include "kadisredu/data_structures/dist_dynamic_graph.h"
#include "kadisredu/data_structures/node_marker.h"
#include "kadisredu/data_structures/sized_vector.h"
#include "kadisredu/definitions.h"
#include "kadisredu/tools/timer.h"

#include <memory>
#include <vector>

namespace kadisredu::mwis {

class Reducer;
using weighted_reduce_ptr = Reducer*;

enum class ReductionType : unsigned {
  neighborhood_removal = 0,
  domination = 1,
  clique = 2,
  simplicial_weight_transfer = 3,
  degree_one = 4,
  degree_two = 5,
  basic_single_edge = 6,
  extended_single_edge = 7,
  generalized_neighborhood_removal = 8,
  neighborhood_folding = 9,
  generalized_neighborhood_folding = 10,
  received_dummy = 11,
  peeler = 12,
  extended_domination = 13
};

static constexpr unsigned MAX_REDUCTION_NUM =
    13;  ////< maximum assigned number among reduction types

inline std::string to_string(const ReductionType& type) {
  switch (type) {
    case ReductionType::neighborhood_removal:
      return "neighborhoood-removal";
    case ReductionType::domination:
      return "domination";
    case ReductionType::clique:
      return "clique";
    case ReductionType::simplicial_weight_transfer:
      return "simplicial-weight-transfer";
    case ReductionType::degree_one:
      return "degree-one";
    case ReductionType::degree_two:
      return "degree-two";
    case ReductionType::basic_single_edge:
      return "basic-single-edge";
    case ReductionType::extended_single_edge:
      return "extended-single-edge";
    case ReductionType::neighborhood_folding:
      return "neighborhood-folding";
    case ReductionType::generalized_neighborhood_folding:
      return "generalized-neighborhood-folding";
    case ReductionType::generalized_neighborhood_removal:
      return "generalized-neighborhood-removal";
    case ReductionType::received_dummy:
      return "received";
    case ReductionType::peeler:
      return "peeler";
    case ReductionType::extended_domination:
      return "extended_domination";
    default:
      return "unknown";
  }
}

class general_weighted_reduction;
using reduction_ptr = std::unique_ptr<general_weighted_reduction>;

class general_weighted_reduction {
 public:
  NodeMarker marker;

  general_weighted_reduction(NodeID num_nodes, NodeID num_ghosts) : marker(num_nodes) {};
  virtual ~general_weighted_reduction() = default;

  virtual ReductionType get_reduction_type() = 0;

  KADISREDU_DEFINE_METHOD_TIMER(reduce_timer);
  virtual bool reduce(weighted_reduce_ptr algo) = 0;

  virtual void restore(weighted_reduce_ptr algo, NodeID modified_node) {};

  virtual void apply(weighted_reduce_ptr algo, NodeID modified_node) {};

  virtual bool reduce_detached_ghost(weighted_reduce_ptr algo, NodeID detached_ghost,
                                     sized_vector<NodeID>& detached_neigh) {
    return false;
  };

  virtual reduction_ptr clone() = 0;
};

template <typename NodeHandler>
inline void for_each_unset_marked_node(const Solution& solution, NodeMarker& marked,
                                       NodeHandler&& f) {
  for (auto node : marked.current()) {
    if (solution.node_status[node] == wis_status::UNSET) {
      f(node);
    }
  }
};

class neighborhood_removal : public general_weighted_reduction {
 public:
  neighborhood_removal(NodeID num_nodes, NodeID num_ghosts)
      : general_weighted_reduction(num_nodes, num_ghosts) {}
  ~neighborhood_removal() override = default;

  ReductionType get_reduction_type() final { return ReductionType::neighborhood_removal; }

  bool reduce(weighted_reduce_ptr algo) override;

  reduction_ptr clone() override { return std::make_unique<neighborhood_removal>(*this); }
};

class generalized_neighborhood_folding : public general_weighted_reduction {
 public:
  generalized_neighborhood_folding(NodeID num_nodes, NodeID num_ghosts)
      : general_weighted_reduction(num_nodes, num_ghosts) {}
  ~generalized_neighborhood_folding() override = default;

  ReductionType get_reduction_type() override {
    return ReductionType::generalized_neighborhood_folding;
  }

  bool reduce(weighted_reduce_ptr algo) override;

  void restore(weighted_reduce_ptr algo, NodeID modified_node) override;
  void apply(weighted_reduce_ptr algo, NodeID modified_node) override;

  reduction_ptr clone() override {
    return std::make_unique<generalized_neighborhood_folding>(*this);
  }

 protected:
  struct fold_nodes {
    NodeID main;
    std::vector<NodeID> MWIS;
  };

  struct restore_data {
    fold_nodes nodes;
    GlobalNodeWeight main_weight;
    GlobalNodeWeight MWIS_weight;
    dist_dynamic_graph::dynamic_neighborhood main_neighbor_list;
    std::vector<std::vector<NodeID>> MWIS_node_vecs;
  };
  std::vector<restore_data> applications;

  void fold(weighted_reduce_ptr algo, NodeID main, fast_set& mwis, GlobalNodeWeight mwis_weight);
};

class neighborhood_folding : public generalized_neighborhood_folding {
 public:
  neighborhood_folding(NodeID num_nodes, NodeID num_ghosts)
      : generalized_neighborhood_folding(num_nodes, num_ghosts) {}
  ~neighborhood_folding() override = default;

  ReductionType get_reduction_type() final { return ReductionType::neighborhood_folding; }

  bool reduce(weighted_reduce_ptr algo) override;
  reduction_ptr clone() override { return std::make_unique<neighborhood_folding>(*this); }
};

class generalized_neighborhood_removal : public general_weighted_reduction {
 public:
  generalized_neighborhood_removal(NodeID num_nodes, NodeID num_ghosts)
      : general_weighted_reduction(num_nodes, num_ghosts) {}
  ~generalized_neighborhood_removal() override = default;

  ReductionType get_reduction_type() final {
    return ReductionType::generalized_neighborhood_removal;
  }

  bool reduce(weighted_reduce_ptr algo) override;

  reduction_ptr clone() override {
    return std::make_unique<generalized_neighborhood_removal>(*this);
  }
};

class degree_one : public general_weighted_reduction {
 public:
  degree_one(NodeID num_nodes, NodeID num_ghosts)
      : general_weighted_reduction(num_nodes, num_ghosts) {}
  ~degree_one() override = default;

  ReductionType get_reduction_type() final { return ReductionType::degree_one; }

  bool reduce(weighted_reduce_ptr algo) override;

  reduction_ptr clone() override { return std::make_unique<degree_one>(*this); }

  void fold(weighted_reduce_ptr algo, NodeID degree_one_node);
  void restore(weighted_reduce_ptr algo, NodeID modified_node) override;
  void apply(weighted_reduce_ptr algo, NodeID modified_node) override;
  bool reduce_detached_ghost(weighted_reduce_ptr algo, NodeID detached_ghost,
                             sized_vector<NodeID>& detached_neigh) override;
  ;

 private:
  struct restore_data {
    NodeWeight offset;
    NodeID node;
  };

  std::vector<restore_data> applications;
};

class extended_domination : public general_weighted_reduction {
 public:
  extended_domination(NodeID num_nodes, NodeID num_ghosts)
      : general_weighted_reduction(num_nodes, num_ghosts) {}
  ~extended_domination() override = default;

  ReductionType get_reduction_type() final { return ReductionType::extended_domination; }

  bool reduce(weighted_reduce_ptr algo) override;

  reduction_ptr clone() override { return std::make_unique<extended_domination>(*this); }

  void fold(weighted_reduce_ptr algo, NodeID v, NodeID u);
  void restore(weighted_reduce_ptr algo, NodeID modified_node) override;
  void apply(weighted_reduce_ptr algo, NodeID modified_node) override;

 private:
  std::vector<NodeID> applications;
};


class v_shape : public general_weighted_reduction {
 public:
  v_shape(NodeID num_nodes, NodeID num_ghosts)
      : general_weighted_reduction(num_nodes, num_ghosts) {}
  ~v_shape() override = default;

  ReductionType get_reduction_type() final { return ReductionType::degree_two; }

  bool reduce(weighted_reduce_ptr algo) override;

  reduction_ptr clone() override { return std::make_unique<v_shape>(*this); }

  void restore(weighted_reduce_ptr algo, NodeID modified_node) override;
  void apply(weighted_reduce_ptr algo, NodeID modified_node) override;

 private:
  struct fold_nodes {
    NodeID main;
    std::array<NodeID, 2> rest;
  };

  struct restore_data {
    fold_nodes nodes;
    NodeWeight main_weight;
    int fold_case;
    dist_dynamic_graph::dynamic_neighborhood main_neighbor_list;
    std::array<std::vector<NodeID>, 2> node_vecs;
  };
  std::vector<restore_data> applications;

  void fold_mid(weighted_reduce_ptr algo, const fold_nodes& nodes);
  void fold_max(weighted_reduce_ptr algo, const fold_nodes& nodes);
};

class domination : public general_weighted_reduction {
 public:
  domination(NodeID num_nodes, NodeID num_ghosts)
      : general_weighted_reduction(num_nodes, num_ghosts) {}
  ~domination() override = default;

  ReductionType get_reduction_type() final { return ReductionType::domination; }

  bool reduce(weighted_reduce_ptr algo) override;

  reduction_ptr clone() override { return std::make_unique<domination>(*this); }
};

class basic_single_edge : public general_weighted_reduction {
 public:
  basic_single_edge(NodeID num_nodes, NodeID num_ghosts)
      : general_weighted_reduction(num_nodes, num_ghosts) {}
  ~basic_single_edge() override = default;

  ReductionType get_reduction_type() final { return ReductionType::basic_single_edge; }

  bool reduce(weighted_reduce_ptr algo) override;

  reduction_ptr clone() override { return std::make_unique<basic_single_edge>(*this); }
};

class extended_single_edge : public general_weighted_reduction {
 public:
  extended_single_edge(NodeID num_nodes, NodeID num_ghosts)
      : general_weighted_reduction(num_nodes, num_ghosts) {}
  ~extended_single_edge() override = default;

  ReductionType get_reduction_type() final { return ReductionType::extended_single_edge; }

  bool reduce(weighted_reduce_ptr algo) override;

  reduction_ptr clone() override { return std::make_unique<extended_single_edge>(*this); }
};

class clique : public general_weighted_reduction {
 public:
  clique(NodeID num_nodes, NodeID num_ghosts) : general_weighted_reduction(num_nodes, num_ghosts) {}
  ~clique() override = default;

  ReductionType get_reduction_type() final { return ReductionType::clique; }

  bool reduce(weighted_reduce_ptr algo) override;

  reduction_ptr clone() override { return std::make_unique<clique>(*this); }
};

class simplicial_weight_transfer : public general_weighted_reduction {
 public:
  simplicial_weight_transfer(NodeID num_nodes, NodeID num_ghosts)
      : general_weighted_reduction(num_nodes, num_ghosts) {}
  ~simplicial_weight_transfer() override = default;

  ReductionType get_reduction_type() final { return ReductionType::simplicial_weight_transfer; }

  bool reduce(weighted_reduce_ptr algo) override;

  reduction_ptr clone() override { return std::make_unique<simplicial_weight_transfer>(*this); }

  void fold(weighted_reduce_ptr algo, NodeID max_simplicial_node);

  void restore(weighted_reduce_ptr algo, NodeID modified_node) override;

  void apply(weighted_reduce_ptr algo, NodeID modified_node) override;

 private:
  struct restore_data {
    NodeWeight offset;
    // sized_vector<NodeID> non_simplicials;
  };

  std::vector<restore_data> applications;
};

template <class Last>
void make_reduction_vector_helper(std::vector<ReductionPtr<general_weighted_reduction>>& vec,
                                  NodeID num_nodes, NodeID num_ghosts) {
  vec.emplace_back(new Last(num_nodes, num_ghosts));
};

template <class First, class Second, class... Redus>
void make_reduction_vector_helper(std::vector<ReductionPtr<general_weighted_reduction>>& vec,
                                  NodeID num_nodes, NodeID num_ghosts) {
  vec.emplace_back(new First(num_nodes, num_ghosts));
  make_reduction_vector_helper<Second, Redus...>(vec, num_nodes, num_ghosts);
};

template <class... Redus>
std::vector<ReductionPtr<general_weighted_reduction>> make_reduction_vector(NodeID num_nodes,
                                                                            NodeID num_ghosts) {
  std::vector<ReductionPtr<general_weighted_reduction>> vec;
  make_reduction_vector_helper<Redus...>(vec, num_nodes, num_ghosts);
  return vec;
};

}  // namespace kadisredu::mwis
