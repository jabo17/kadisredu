//
// Created by jannickb on 6/17/24.
//

#pragma once

#include "kadisredu/data_structures/fast_set.h"
#include "kadisredu/data_structures/sized_vector.h"
#include "kadisredu/definitions.h"
#include "kadisredu/tools/logger.h"
#include "kadisredu/tools/utils.h"

#include <absl/container/flat_hash_map.h>
#include <algorithm>
#include <google/dense_hash_map>
#include <iterator>
#include <kagen/kagen.h>
#include <kaminpar.h>
#include <kamping/checking_casts.hpp>
#include <kamping/collectives/allgather.hpp>
#include <kamping/collectives/allreduce.hpp>
#include <kamping/collectives/alltoall.hpp>
#include <kamping/communicator.hpp>
#include <kamping/data_buffer.hpp>
#include <kamping/environment.hpp>
#include <kamping/named_parameters.hpp>
#include <kamping/utils/flatten.hpp>
#include <stack>

namespace kadisredu {

/***
 * @class dist_dynamic_graph represents the local block of a distributed simple,
 *  undirected, node-weighted graph while allowing to remove nodes, edges
 *  and restore them in the reverse order.
 * @brief: The graph is limited as weights need to be representable by
 * kadisredu::GlobalNodeWeight, node ids (without ghosts) by kadisredu::NodeID,
 * and global node ids by kadisredu::GlobalNodeID, respectively.
 */
class dist_dynamic_graph {
 public:
  using BlockID = kadisredu::BlockID;
  using NodeID = kadisredu::NodeID;
  using EdgeID = kadisredu::EdgeID;
  using NodeWeight = kadisredu::NodeWeight;
  using GlobalNodeID = kadisredu::GlobalNodeID;
  using GlobalEdgeID = kadisredu::GlobalEdgeID;
  using GlobalNodeWeight = kadisredu::GlobalNodeWeight;

  class dynamic_neighborhood {
    friend dist_dynamic_graph;

   public:
    dynamic_neighborhood() = default;
    explicit dynamic_neighborhood(std::span<NodeID> neighbors) { build(neighbors); }
    void build(std::vector<NodeID>&& neighbors) {
      neigh = std::move(neighbors);
      visible = neigh.size();
    }
    void build(std::span<NodeID> neighbors) {
      neigh = std::vector<NodeID>(neighbors.begin(), neighbors.end());
      visible = neigh.size();
    }

    [[nodiscard]] NodeID deg() const { return visible; }
    [[nodiscard]] NodeID size() const { return visible; }

    NodeID hidden() { return neigh.size() - visible; }

    using iterator = typename std::vector<NodeID>::iterator;
    using const_iterator = typename std::vector<NodeID>::const_iterator;
    using reverse_iterator = typename std::reverse_iterator<iterator>;
    [[nodiscard]] constexpr iterator begin() noexcept { return neigh.begin(); }
    [[nodiscard]] constexpr const_iterator begin() const noexcept { return neigh.begin(); }
    [[nodiscard]] constexpr const_iterator cbegin() const noexcept { return neigh.cbegin(); }
    [[nodiscard]] constexpr iterator hidden_end() noexcept { return neigh.end(); }
    [[nodiscard]] constexpr const_iterator hidden_end() const noexcept { return neigh.end(); }
    [[nodiscard]] constexpr iterator end() noexcept { return neigh.begin() + visible; }
    [[nodiscard]] constexpr const_iterator end() const noexcept { return neigh.begin() + visible; }
    [[nodiscard]] constexpr const_iterator cend() const noexcept {
      return neigh.cbegin() + visible;
    }

    [[nodiscard]] constexpr reverse_iterator rbegin() noexcept {
      return std::make_reverse_iterator(end());
    }
    [[nodiscard]] constexpr reverse_iterator rend() noexcept {
      return std::make_reverse_iterator(begin());
    }

    [[nodiscard]] NodeID& operator[](EdgeID edge) { return neigh[edge]; }

    [[nodiscard]] const NodeID& operator[](EdgeID edge) const { return neigh[edge]; }

    [[nodiscard]] NodeID& at(EdgeID edge) {
      // TODO range check
      return neigh[edge];
    }

    [[nodiscard]] const NodeID& at(const EdgeID edge) const { return neigh.at(edge); }

    [[nodiscard]] NodeID ghost_deg() const { return ghosts; }

   private:
    void restore_last() {
      KASSERT(visible < neigh.size());
      KASSERT(visible == std::distance(begin(), end()));

      ++visible;

      KASSERT(visible == std::distance(begin(), end()));
    }

    void hide(iterator iter) {
      KASSERT(iter >= begin());
      KASSERT(iter < end());
      KASSERT(visible > 0);

      using std::swap;
      swap(*iter, neigh[--visible]);

      KASSERT(visible >= 0);
      KASSERT(visible == std::distance(begin(), end()));
    }

    std::vector<NodeID> neigh;
    NodeID visible = 0;
    NodeID ghosts = 0;
  };

  class internal_ghost_node_map_old {
   public:
    using Map = google::dense_hash_map<GlobalNodeID, NodeID>;

    std::vector<BlockID> blocks;
    std::vector<GlobalNodeID> global_ids;
    Map local_id_map;
    static_assert(std::is_move_constructible<Map>::value, "Map must be move constructible");

    explicit internal_ghost_node_map_old(NodeID reserve = 0) : local_id_map(reserve) {
      blocks.reserve(reserve);
      global_ids.reserve(reserve);

      local_id_map.set_deleted_key(std::numeric_limits<GlobalNodeID>::max());
      local_id_map.set_empty_key(std::numeric_limits<GlobalNodeID>::max() - 1);
      local_id_map.min_load_factor(0.0);
    };

    std::pair<NodeID, bool> add(GlobalNodeID global_id, BlockID block_id) {
      auto [it, succ] = local_id_map.insert({global_id, global_ids.size()});
      KASSERT(local_id_map.find(global_id) != local_id_map.end());
      if (succ) {
        global_ids.push_back(global_id);
        blocks.push_back(block_id);
        KASSERT(get(global_id) == global_ids.size() - 1);
      }
      return {it->second, succ};
    }

    [[nodiscard]] NodeID get(GlobalNodeID global_id) {
      KASSERT(local_id_map.find(global_id) != local_id_map.end());
      return local_id_map[global_id];
    }
    [[nodiscard]] NodeID num_ghosts() const {
      KASSERT(blocks.size() == global_ids.size());
      KASSERT(local_id_map.size() == global_ids.size());
      return global_ids.size();
    }
  };

  class internal_ghost_node_map {
   public:
    using Map = absl::flat_hash_map<GlobalNodeID, NodeID>;

    std::vector<BlockID> blocks;
    std::vector<GlobalNodeID> global_ids;
    Map local_id_map;
    static_assert(std::is_move_constructible<Map>::value, "Map must be move constructible");

    explicit internal_ghost_node_map(NodeID reserve = 0) {
      local_id_map.reserve(reserve);
      blocks.reserve(reserve);
      global_ids.reserve(reserve);
    };

    std::pair<NodeID, bool> add(GlobalNodeID global_id, BlockID block_id) {
      auto [it, succ] = local_id_map.insert({global_id, global_ids.size()});
      KASSERT(local_id_map.find(global_id) != local_id_map.end());
      if (succ) {
        global_ids.push_back(global_id);
        blocks.push_back(block_id);
        KASSERT(get(global_id) == global_ids.size() - 1);
      }
      return {it->second, succ};
    }

    [[nodiscard]] NodeID get(GlobalNodeID global_id) {
      KASSERT(local_id_map.find(global_id) != local_id_map.end());
      return local_id_map.find(global_id)->second;
    }

    [[nodiscard]] BlockID get_block_from_global(GlobalNodeID global_id) {
      return blocks[get(global_id)];
    }

    [[nodiscard]] NodeID num_ghosts() const {
      KASSERT(blocks.size() == global_ids.size());
      KASSERT(local_id_map.size() == global_ids.size());
      return global_ids.size();
    }
  };

  template <typename Update>
  struct ghost_update {
    GlobalNodeID global_id;
    Update update;
  };

  using ghost_node_mapping = internal_ghost_node_map;

  /**
   * @brief Construct a Distributed Dynamic Graph.
   * @param comm MPI Communicator
   * @param first_global_node first global node ID in distributed graph
   * @param number_nodes amount of local nodes
   * @param ghosts_mapping ghost node mapper mapping global node IDs of ghosts to local node IDs to
   * [num_nodes, num_nodes+|ghosts|)
   */
  dist_dynamic_graph(const MPI_Comm comm, const GlobalNodeID first_global_node,
                     const NodeID number_nodes,
                     dist_dynamic_graph::ghost_node_mapping ghosts_mapping)
      : comm(comm),
        first_global_node(first_global_node),
        number_nodes(number_nodes),
        visible_nodes(number_nodes),
        visible_ghosts(ghosts_mapping.num_ghosts()),
        border(number_nodes),
        border_set(number_nodes),
        ghosts(std::move(ghosts_mapping)),
        adjs(number_nodes + ghosts.num_ghosts()),
        hidden(number_nodes + ghosts.num_ghosts()),
        hidden_nodes(number_nodes + ghosts.num_ghosts()),
        weights(number_nodes + ghosts.num_ghosts()),
        bulk_buffer(number_nodes + ghosts.num_ghosts()),
        bulk_buffer_set(number_nodes + ghosts.num_ghosts()) {
    KASSERT(num_nodes() + num_ghosts() == adjs.size());
  }

  dist_dynamic_graph() : number_nodes(0), border_set(0), bulk_buffer_set(0) {};

  dist_dynamic_graph(const dist_dynamic_graph& other, std::vector<NodeID>& mapping)
      : comm(other.comm),
        number_nodes(other.visible_nodes),
        visible_nodes(other.visible_nodes),
        visible_ghosts(other.visible_ghosts),
        adjs(other.visible_nodes + other.visible_ghosts),
        hidden(other.visible_nodes + other.visible_ghosts),
        hidden_nodes(other.visible_nodes + other.visible_ghosts),
        weights(other.visible_nodes + other.visible_ghosts),
        border_set(other.visible_nodes),
        border(other.visible_nodes),
        bulk_buffer(other.visible_nodes + other.visible_ghosts),
        bulk_buffer_set(other.visible_nodes + other.visible_ghosts) {
    copy_visible(other, mapping);
  }

  /**
   *  @return number of blocks (number of processes which store a part of this graph)
   */
  BlockID blocks() { return comm.size(); }

  /**
   * @return amount of initial nodes at this rank excluding ghosts
   */
  [[nodiscard]] NodeID num_nodes() const { return number_nodes; }

  /**
   * @return amount of visible nodes at this rank excluding ghosts
   */
  [[nodiscard]] NodeID num_visible_nodes() const { return visible_nodes; }

  [[nodiscard]] bool is_visible(NodeID node) const {
    KASSERT(node < num_nodes() + num_ghosts());
    return !hidden[node];
  }

  /**
   * @return amount of initial ghosts
   */
  [[nodiscard]] NodeID num_ghosts() const { return ghosts.num_ghosts(); }

  [[nodiscard]] NodeID num_visible_ghosts() const { return visible_ghosts; }

  /**
   * @param call function that is called on a node id of a visible non-ghost
   * node
   */
  template <typename UnaryFunc>
  void for_each_visible_node(UnaryFunc&& func) const {
    for (NodeID node = 0; node < num_nodes(); ++node) {
      if (!hidden[node]) {
        func(node);
      }
    }
  };

  template <typename UnaryFunc>
  void for_each_visible_node_inc_ghosts(UnaryFunc&& func) const {
    for (NodeID node = 0; node < num_nodes() + num_ghosts(); ++node) {
      if (!hidden[node]) {
        func(node);
      }
    }
  };

  template <typename UnaryFunc>
  void for_each_visible_ghost(UnaryFunc&& func) const {
    for (NodeID node = num_nodes(); node < num_nodes() + num_ghosts(); ++node) {
      if (!hidden[node]) {
        func(node);
      }
    }
  };

  /**
   * @brief hide a visible node
   * @brief hides all incoming edges to this node
   * @param iter iterator that points to a node in visible_nodes()
   */
  void hide_node(NodeID node) {
    KASSERT(node < num_nodes() + num_ghosts());
    KASSERT(node >= 0);
    KASSERT(!hidden[node]);

    // set hidden
    hidden[node] = true;
    hidden_nodes.push_back(node);
    visible_nodes -= static_cast<int>(node < num_nodes());
    visible_ghosts -= static_cast<int>(node >= num_nodes());
    // hide incoming edges
    auto& node_neigh = adjs[node];
    bool node_is_ghost = node >= num_nodes();
    std::for_each(node_neigh.begin(), node_neigh.end(), [node, node_is_ghost, this](NodeID target) {
      KASSERT(target != node);
      KASSERT(!hidden[target]);
      auto& target_neigh = this->adjs[target];
      // find incoming edges by scanning over neighborhoods
      KASSERT(std::find(target_neigh.begin(), target_neigh.end(), node) != target_neigh.end());
      target_neigh.hide(std::find(target_neigh.begin(), target_neigh.end(), node));
      if (node_is_ghost) {
        target_neigh.ghosts--;
        KASSERT(target_neigh.ghosts ==
                std::count_if(target_neigh.begin(), target_neigh.end(),
                              [&](NodeID target2) { return is_ghost(target2); }));
      }
    });

    KASSERT(hidden_nodes.back() == node);

    KASSERT(hidden[node]);
  }

  void hide_edge(NodeID node, NodeID target) {
    auto& neigh = this->adjs[node];
    // find incoming edges by scanning over neighborhoods
    KASSERT(std::find(neigh.begin(), neigh.end(), target) != neigh.end());
    neigh.hide(std::find(neigh.begin(), neigh.end(), target));
    if (is_ghost(target)) {
      neigh.ghosts--;
      KASSERT(neigh.ghosts == std::count_if(neigh.begin(), neigh.end(),
                                            [&](NodeID target2) { return is_ghost(target2); }));
    }
  }

  void bulk_hide(std::span<const NodeID> nodes, const fast_set& nodes_set) {
    bulk_buffer.clear();
    bulk_buffer_set.clear();
    NodeID new_hidden_ghosts = 0;
    for (auto& node : nodes) {
      hidden[node] = true;
      hidden_nodes.push_back(node);
      visible_nodes -= static_cast<int>(node < num_nodes());
      visible_ghosts -= static_cast<int>(node >= num_nodes());

      auto& node_neigh = adjs[node];

      // hide edges to hidden vertices of the bulk
      auto new_end = utils::partition(node_neigh, [&](NodeID target) { return hidden[target]; });
      node_neigh.visible = std::distance(node_neigh.begin(), new_end);
      node_neigh.ghosts -= new_hidden_ghosts;
      KASSERT(node_neigh.ghosts == std::count_if(node_neigh.begin(), node_neigh.end(),
                                                 [&](NodeID target) { return is_ghost(target); }));

      new_hidden_ghosts += static_cast<int>(node >= num_nodes());

      std::for_each(node_neigh.begin(), node_neigh.end(), [&](NodeID target) {
        KASSERT(!hidden[target]);
        if (!nodes_set.contains(target) && bulk_buffer_set.add(target)) {
          bulk_buffer.push_back(target);
        }
      });
    }

    for (auto& node : bulk_buffer) {
      auto& node_neigh = adjs[node];
      // bulk hide all edges from modified nodes to the hidden nodes
      auto new_end = utils::partition(node_neigh, [&](NodeID target) { return hidden[target]; });
      node_neigh.ghosts -= std::count_if(new_end, node_neigh.end(),
                                         [&](NodeID target) { return target >= num_nodes(); });
      node_neigh.visible = std::distance(node_neigh.begin(), new_end);
      KASSERT(node_neigh.ghosts == std::count_if(node_neigh.begin(), node_neigh.end(),
                                                 [&](NodeID target) { return is_ghost(target); }));
    }
  }

  /**
   * @brief restore the last hidden node and incoming edges
   */
  NodeID restore_last() {
    KASSERT(num_visible_nodes() >= 0);
    KASSERT(num_visible_nodes() < num_nodes() + num_ghosts());

    NodeID restored_node = hidden_nodes.back();
    hidden_nodes.pop_back();
    KASSERT(hidden[restored_node] == true);
    hidden[restored_node] = false;
    visible_nodes += static_cast<NodeID>(restored_node < num_nodes());
    visible_ghosts += static_cast<NodeID>(restored_node >= num_nodes());
    // restore incoming edges
    bool node_is_ghost = restored_node >= num_nodes();
    std::for_each(adjs[restored_node].begin(), adjs[restored_node].end(),
                  [node_is_ghost, this](NodeID target) {
                    auto& neigh = this->adjs[target];
                    neigh.restore_last();
                    if (node_is_ghost) {
                      neigh.ghosts++;
                    }
                  });
    return restored_node;
  }

  void restore_and_replace_edge(NodeID node, NodeID new_neighbor) {
    auto& slot = adjs[node];
    KASSERT(slot.visible < slot.neigh.size());
    slot.neigh[slot.visible++] = new_neighbor;
  }

  void replace_last_hidden_edge(NodeID node, NodeID new_neighbor) {
    auto& slot = adjs[node];
    KASSERT(slot.visible < slot.neigh.size());
    slot.neigh[slot.visible] = new_neighbor;
  }

  void restore_last_hidden_edge(NodeID node) {
    auto& slot = adjs[node];
    KASSERT(slot.visible < slot.neigh.size());
    slot.restore_last();
    if (is_ghost(slot[slot.visible - 1])) {
      slot.ghosts++;
    }
  }

  void add_edge(NodeID node, NodeID new_neighbor) {
    // based on m^2wis+s by Grossman et al
    auto& slot = adjs[node];
    KASSERT(std::find(slot.begin(), slot.end(), new_neighbor) == slot.end());
    slot.neigh.push_back(new_neighbor);
    slot.visible++;
    for (std::size_t i = slot.neigh.size() - 1; i >= slot.visible; i--) {
      using std::swap;
      swap(slot.neigh[i - 1], slot.neigh[i]);
    }
  }

  void remove_hidden_edge(NodeID node, NodeID neighbor) {
    // based on m^2wis+s by Grossman et al
    auto& slot = adjs[node];
    KASSERT(std::find(slot.end(), slot.hidden_end(), neighbor) != slot.hidden_end());
    for (std::size_t i = slot.visible; i < slot.neigh.size() - 1; i++) {
      if (slot.neigh[i] == neighbor) {
        using std::swap;
        swap(slot.neigh[i], slot.neigh[i + 1]);
      }
    }
    slot.neigh.pop_back();
  }

  /**
   * @brief get the weight of a node
   * @brief node might be hidden
   * @param node some node id of [0, num_nodes()+num_ghosts())
   */
  [[nodiscard]] NodeWeight get_weight(NodeID node) const { return weights[node]; }
  void set_weight(NodeID node, NodeWeight weight) { weights[node] = weight; }

  [[nodiscard]] NodeWeight get_weight_neigh(NodeID node) const {
    KASSERT(!is_ghost(node));
    return std::transform_reduce(
        adjs[node].begin(), adjs[node].end(), (NodeWeight)0,
        [](NodeWeight partial_sum, NodeWeight weight) { return partial_sum + weight; },
        [this](NodeID target) { return get_weight(target); });
  }

  [[nodiscard]] NodeWeight get_max_weight_neigh(NodeID node) const {
    KASSERT(!is_ghost(node));
    return std::transform_reduce(
        adjs[node].begin(), adjs[node].end(), (NodeWeight)0,
        [](NodeWeight current_max, NodeWeight weight) { return std::max(current_max, weight); },
        [this](NodeID target) { return get_weight(target); });
  }

  /**
   * @brief obtain const reference to a dynamic neighborhood
   * @brief dynamic neighborhood hides all neighbors that were hided before node
   *   was
   * @param node some node id of [0, num_nodes())
   */
  const dynamic_neighborhood& operator[](NodeID node) const {
    KASSERT(node < num_nodes() + num_ghosts());
    return adjs[node];
  }

  /**
   * @brief build a neighborhood for some node
   * @brief assumes ghost mapper ist set
   */
  void build_neighborhood(NodeID node, std::vector<NodeID>&& neighbors) {
    KASSERT(node < adjs.size());
    KASSERT(neighbors.size() == 0 ||
            *std::max_element(neighbors.begin(), neighbors.end()) < num_nodes() + num_ghosts());
    adjs[node].build(std::move(neighbors));
  }
  void build_neighborhood(NodeID node, std::span<NodeID> neighbors) {
    KASSERT(node < adjs.size());
    KASSERT(neighbors.size() == 0 ||
            *std::max_element(neighbors.begin(), neighbors.end()) < num_nodes() + num_ghosts());
    adjs[node].build(neighbors);
  }
  void set_neighborhood(NodeID node, dist_dynamic_graph::dynamic_neighborhood&& neighbors) {
    KASSERT(node < adjs.size());
    KASSERT(neighbors.neigh.size() == 0 ||
            *std::max_element(neighbors.neigh.begin(), neighbors.neigh.end()) <
                num_nodes() + num_ghosts());
    adjs[node] = std::move(neighbors);
  }

  /**
   * @brief test whether some node is a ghost
   * @param node some node id of [0, num_nodes()+num_ghosts())
   * @return true iff. node is not less than num_nodes()
   */
  [[nodiscard]] bool is_ghost(NodeID node) const {
    KASSERT(node < num_nodes() + num_ghosts());
    return node >= num_nodes();
  }

  /**
   * @brief given a ghost by its global node ID get the node ID to which it is
   *   mapped
   * @param global_id global ID of ghost in distributed graph
   * @return ghost node ID at this rank in [num_nodes(),
   *   num_nodes()+num_ghosts())
   */
  [[nodiscard]] NodeID get_ghost(GlobalNodeID global_id) const {
    KASSERT(ghosts.local_id_map.find(global_id) != ghosts.local_id_map.end());
    return (*ghosts.local_id_map.find(global_id)).second + num_nodes();
  }

  [[nodiscard]] BlockID get_ghost_block(NodeID node) const {
    KASSERT(is_ghost(node));
    return ghosts.blocks[node - num_nodes()];
  }

  /**
   * @return first global node id that is represented in this local block
   */
  [[nodiscard]] GlobalNodeID get_first_global_node() const { return first_global_node; }

  /**
   * @brief get the global node id
   * @param node given a node of range [0,num_nodes()+num_ghosts())
   * @return the global node id in the distributed graph
   */
  [[nodiscard]] GlobalNodeID get_global_node_id(NodeID node) const {
    KASSERT(node >= 0);
    KASSERT(node < num_nodes() + num_ghosts());
    if (node < num_nodes()) {
      return get_first_global_node() + node;
    } else {
      return ghosts.global_ids[node - num_nodes()];
    }
  }

  /**
   * @brief Helper to build sequence of border vertices and ghost neighborhoods
   */
  void build_border() {
    border.clear();
    border_set.clear();
    if (num_ghosts() == 0) {
      return;
    }
    // build sequence of border vertices and count their ghost neighbors
    for (NodeID node = 0; node < num_nodes(); ++node) {
      for (auto& target : adjs[node]) {
        if (is_ghost(target)) {
          if (border_set.add(node)) {
            border.push_back(node);
          }
          KASSERT(target - num_nodes() < num_ghosts());
          adjs[node].ghosts++;
          adjs[target].visible++;
        }
      }
      KASSERT(adjs[node].ghosts == std::count_if(adjs[node].begin(), adjs[node].end(),
                                                 [&](NodeID target) { return is_ghost(target); }));
    }
    // allocate memory for ghost neighborhoods
    for (NodeID ghost = 0; ghost < num_ghosts(); ++ghost) {
      auto& slot = adjs[ghost + num_nodes()];
      slot.neigh.resize(slot.visible);
      slot.visible = 0;
      KASSERT(adjs[ghost + num_nodes()].neigh.size() >= 1,
              "Each partial ghost should be adjacent to at least one border node.");
    }
    // build ghost neighborhoods
    for (NodeID node : border) {
      for (auto& target : adjs[node]) {
        if (is_ghost(target)) {
          auto& slot = adjs[target];
          slot.neigh[slot.visible++] = node;
        }
      }
    }
  }

  /**
   * @brief obtain node id of nodes that were originally (before hiding) at the
   *   border (adjacent to a ghost)
   * @return const reference to a container holding these node ids
   */
  [[nodiscard]] const sized_vector<NodeID>& get_border() const { return border; }

  /**
   * @brief check whether a node is originally (before hiding) a border node
   * @param node the node id which should be checked
   * @return true iff. it is a border node
   */
  [[nodiscard]] bool is_original_border_node(NodeID node) const {
    KASSERT(node < num_nodes() + num_ghosts());
    if (node < num_nodes()) {
      return border_set.contains(node);
    } else {
      return false;
    }
  }

  [[nodiscard]] bool is_active_border_node(NodeID node) const {
    KASSERT(node < num_nodes() + num_ghosts());
    if (node < num_nodes()) {
      return adjs[node].ghost_deg() > 0;
    } else {
      return false;
    }
  }

  /**
   * @brief copy the visible part of a graph (forget hidden parts)
   * @brief assumes graph is in a state where interfaces match
   * @param other graph that is copied
   * @param mapping maps new node ids to old ones; asserts that the sized_vector
   *   has capacity of at least other.num_visible_nodes() +
   * other.num_visible_ghosts()
   */
  void copy_visible(const dist_dynamic_graph& other, std::vector<NodeID>& mapping) {
    using namespace kamping;
#ifndef NDEBUG
    // check ghost interface
    // all visible border nodes are sent to their adjacent visible ghost partitions
    // we expect that all received ghosts are visible border nodes at each partition
    // sized_vector<std::pair<NodeID, kabool>> updates(
    //     std::min(other.get_border().size(),
    //     static_cast<std::size_t>(other.num_visible_nodes())));
    // std::for_each(other.get_border().begin(), other.get_border().end(),
    //               [&other, &updates](NodeID border_node) {
    //                 if (other.is_visible(border_node)) {
    //                   updates.push_back({border_node, true});
    //                 }
    //               });
    // auto recv_ghosts = other.sync_ghosts(updates);
    // for (auto recv_ghost : recv_ghosts) {
    //   KASSERT(other.is_visible(other.get_ghost(recv_ghost.global_id)), "ghost is hidden.");
    // }
    // KASSERT(recv_ghosts.size() == other.num_visible_ghosts(),
    //         "amount of received visible ghosts does not match local visible ghosts.");
#endif
    comm = other.comm;

    std::vector<NodeID> reverse_mapping(other.num_nodes() + other.num_ghosts());
    mapping.clear();
    // contruct mapping between other's nodes and visible nodes
    other.for_each_visible_node_inc_ghosts([&mapping, &reverse_mapping](auto visible_node) {
      reverse_mapping[visible_node] = mapping.size();
      mapping.push_back(visible_node);
    });
    KASSERT(mapping.size() == other.num_visible_nodes() + other.num_visible_ghosts());

    std::vector<GlobalNodeID> visible_vtx_counts(comm.size(), 0);
    visible_vtx_counts[comm.rank()] = other.num_visible_nodes();
    comm.allgather_inplace(send_recv_buf(visible_vtx_counts));

    first_global_node = std::reduce(visible_vtx_counts.begin(),
                                    visible_vtx_counts.begin() + comm.rank(), (GlobalNodeID)0);
    number_nodes = other.num_visible_nodes();
    visible_nodes = number_nodes;
    visible_ghosts = other.num_visible_ghosts();

    // get new global ids of visible ghosts
    sized_vector<std::pair<NodeID, GlobalNodeID>> border_updates(
        other.border.size());  // upper bound
    std::for_each(
        other.border.begin(), other.border.end(),
        [&other, &border_updates, &reverse_mapping, this](NodeID node) {
          if (other.is_visible(node)) {
            return border_updates.push_back({node, reverse_mapping[node] + first_global_node});
          }
        });
    auto new_global_ghost_ids = other.sync_ghosts(border_updates);
    KASSERT(visible_ghosts == new_global_ghost_ids.size(),
            "Maybe the interfaces do not match (possibly due to isolated ghosts).");
    ghosts.global_ids.resize(visible_ghosts);
    ghosts.blocks.resize(visible_ghosts);
    for (auto [old_global_node_id, new_global_node_id] : new_global_ghost_ids) {
      auto old_ghost = other.get_ghost(old_global_node_id);
      KASSERT(other.is_visible(old_ghost));
      // add ghosts
      // NOTE: we cannot simply use ghosts.add because ghost updates might not be sorted by
      // old_global_node_id
      auto ghost_entry_pos = reverse_mapping[old_ghost] - num_nodes();
      ghosts.local_id_map.insert({new_global_node_id, ghost_entry_pos});
      ghosts.global_ids[ghost_entry_pos] = new_global_node_id;
      ghosts.blocks[ghost_entry_pos] = other.ghosts.blocks[old_ghost - other.num_nodes()];
    }
    KASSERT(visible_ghosts == num_ghosts());

    // resizing
#if KASSERT_ENABLED(KASSERTION_LEVEL_NORMAL)
    std::ranges::for_each(hidden, [](auto h) { KASSERT(!h); });
#endif
    hidden.resize(number_nodes + visible_ghosts, false);
    hidden_nodes.resize(number_nodes + visible_ghosts);
    hidden_nodes.clear();
    border_set.resize(number_nodes);
    weights.resize(number_nodes + visible_ghosts);
    adjs.resize(number_nodes + visible_ghosts);
    KASSERT(adjs.size() == num_nodes() + num_ghosts());

    // build neighborhoods
    NodeID current = 0;
    for (auto node : mapping) {
      KASSERT(other.is_visible(node));
      std::vector<NodeID> new_neighborhood;
      new_neighborhood.reserve(other[node].visible);
      std::ranges::transform(other[node], std::back_inserter(new_neighborhood),
                             [&reverse_mapping](NodeID target) { return reverse_mapping[target]; });
      if (node >= other.num_nodes()) {
        // node is ghost
        KASSERT(other.is_ghost(node));
        // neighbors are border vertices
        for (auto border_node : new_neighborhood) {
          if (border_set.add(border_node)) {
            border.push_back(border_node);
          }
          adjs[border_node].ghosts++;
        }
      }
      build_neighborhood(current, std::move(new_neighborhood));
      ++current;
    }
    KASSERT(current == number_nodes + visible_ghosts);

    // build mapping
    current = 0;
    for (auto node : mapping) {
      weights[current++] = other.weights[node];
    }
  }

  /**
   * @brief build a static reduced graph in csr representation.
   */
  [[nodiscard]] std::pair<kagen::Graph, std::vector<NodeID>> build_visible_graph() const {
    using namespace kagen;
    using namespace kamping;

    Graph graph;
    std::vector<NodeID> reverse_mapping(num_nodes() + num_ghosts(), 0);
    graph.representation = GraphRepresentation::CSR;
    auto& xadj = graph.xadj;  // i-th field stores first index to first edges in adjncy, i+1-the
                              // stores the index of the last edge + 1
    auto& adjncy = graph.adjncy;
    auto& vwgt = graph.vertex_weights;
    auto& adjwgt = graph.edge_weights;

    xadj.resize(num_visible_nodes() + 1);
    vwgt.resize(num_visible_nodes());

    SInt edges = 0;
    SInt nodes_inc_ghosts = 0;
    for_each_visible_node([&edges, &reverse_mapping, &nodes_inc_ghosts, this](NodeID node) {
      edges += adjs[node].deg();
      reverse_mapping[node] = nodes_inc_ghosts++;
    });
    for_each_visible_ghost([&edges, &reverse_mapping, &nodes_inc_ghosts, this](NodeID node) {
      reverse_mapping[node] = nodes_inc_ghosts++;
    });
    adjncy.resize(edges);

    std::vector<GlobalNodeID> node_counts(comm.size(), 0);
    std::vector<GlobalNodeID> node_dist(comm.size(), 0);
    node_counts[comm.rank()] = visible_nodes;
    comm.allgather_inplace(send_recv_buf(node_counts));
    std::inclusive_scan(node_counts.begin(), node_counts.end(), node_dist.begin());
    GlobalNodeID reduced_first_global_node = node_dist[comm.rank()] - visible_nodes;
    graph.vertex_range = std::make_pair(reduced_first_global_node, node_dist[comm.rank()]);

    sized_vector<std::pair<NodeID, GlobalNodeID>> border_updates(border.size());  // upper bound
    std::for_each(
        border.begin(), border.end(),
        [&border_updates, &reverse_mapping, &reduced_first_global_node, this](NodeID node) {
          if (is_visible(node)) {
            KASSERT(reverse_mapping[node] < visible_nodes);
            return border_updates.push_back(
                {node, reverse_mapping[node] + reduced_first_global_node});
          }
        });
    auto new_global_ghost_ids = sync_ghosts(border_updates);
    KASSERT(visible_ghosts == new_global_ghost_ids.size(),
            "Maybe the interfaces do not match (possibly due to isolated ghosts).");

    ghost_node_mapping reduced_ghosts(visible_ghosts);
    for (auto [old_global_node_id, new_global_node_id] : new_global_ghost_ids) {
      auto old_ghost = get_ghost(old_global_node_id);
      KASSERT(is_visible(old_ghost));
      KASSERT(adjs[old_ghost].deg() > 0);
      // add ghosts
      // note: we cannot simply use 'add' of the ghost mapper
      auto ghost_entry_pos = reverse_mapping[old_ghost] - visible_nodes;
      // reduced_ghosts.local_id_map.insert({new_global_node_id, ghost_entry_pos});
      reduced_ghosts.global_ids[ghost_entry_pos] = new_global_node_id;
      reduced_ghosts.blocks[ghost_entry_pos] = ghosts.blocks[old_ghost - num_nodes()];
    }

    SInt new_node = 0;
    xadj[0] = 0;
    edges = 0;
    for_each_visible_node([&xadj, &adjncy, &vwgt, &adjwgt, &new_node, &reduced_ghosts,
                           &reverse_mapping, &reduced_first_global_node, &edges, &node_dist,
                           this](NodeID node) {
      KASSERT(new_node + 1 < xadj.size());
      for (auto target : adjs[node]) {
        KASSERT(edges < adjncy.size());
        SInt& global_target = adjncy[edges];
        if (is_ghost(target)) {
          global_target = reduced_ghosts.global_ids[reverse_mapping[target] - visible_nodes];
          KASSERT(global_target < node_dist[comm.size() - 1]);
        } else {
          global_target = reverse_mapping[target] + reduced_first_global_node;
          KASSERT(global_target < node_dist[comm.rank()]);
        }
        edges++;
      }
      vwgt[new_node++] = get_weight(node);
      xadj[new_node] = edges;
    });

    return std::make_pair(std::move(graph), std::move(reverse_mapping));
  }

  /**
   * @brief send updates concerning border vertices to the respective blocks of
   *   visible ghosts using an alltoallv.
   * @details The update of a border vertex b at block A (=adj. to
   *   other blocks != A) reaches all blocks of ghosts that were adjacent to b until
   *   removing b.
   *   Consequently, if a border vertex b2, located at another block, does not
   *   receive an update from b, it means either:
   *      - there was no update concerning b or
   *      - block A assumes b2 is hidden when the update concerning b was written.
   * @param border_updates contains the updates where
   *   the first value is the local border node and the second value is the
   *   update of type Update
   * @tparam Update a type holding the update information for a single border
   *   node
   * @return a vector of update for the ghost vertices of type
   *   ghost_update<Upate>; the ghost can be identified by the global_id (local
   *   ids can be obtained with get_ghost)
   */
  template <typename Update>
  std::vector<ghost_update<Update>> sync_ghosts(
      sized_vector<std::pair<NodeID, Update>>& border_updates) const {
    return sync_ghosts(border_updates, [](NodeID node) { return true; });
  }
  template <typename Update, typename UnaryFunc>
  std::vector<ghost_update<Update>> sync_ghosts(
      sized_vector<std::pair<NodeID, Update>>& border_updates, UnaryFunc filter_updates) const {
    using namespace kamping;

    std::vector<std::vector<ghost_update<Update>>> buckets(comm.size());
    fast_set reached_blocks_by_border_node(comm.size());

    // scan updated border nodes and prepare update messages for each reached
    // blocks
    for (auto& b_up : border_updates) {
      if (!filter_updates(b_up.first)) {
        continue;
      }
      if (is_ghost(b_up.first)) {
        auto ghost = b_up.first;
        buckets[ghosts.blocks[ghost - num_nodes()]].push_back(
            {.global_id = ghosts.global_ids[ghost - num_nodes()], .update = b_up.second});
      } else {
        reached_blocks_by_border_node.clear();
        auto border_node = b_up.first;
        KASSERT(is_original_border_node(border_node));
        for (auto target : adjs[border_node]) {
          if (is_ghost(target)) {
            auto block = ghosts.blocks[target - num_nodes()];
            KASSERT(block < comm.size());
            if (reached_blocks_by_border_node.add(block)) {
              buckets[block].push_back(
                  {.global_id = border_node + get_first_global_node(), .update = b_up.second});
            }
          }
        }
      }
    }

    // return received ghost updates
    return with_flattened(buckets, comm.size()).call([&](auto... flattened) {
      return comm.alltoallv(std::move(flattened)...);
    });
  };

  template <typename Update, typename UnaryFunc>
  std::vector<ghost_update<Update>> sync_ghosts_with_lazy_updates(
      std::span<const NodeID> update_nodes, UnaryFunc&& compute_update) const {
    using namespace kamping;

    std::vector<std::vector<ghost_update<Update>>> buckets(comm.size());
    fast_set reached_blocks_by_border_node(comm.size());

    // scan updated border nodes and prepare update messages for each reached
    // blocks
    for (const auto& v : update_nodes) {
      std::optional<Update> update = compute_update(v);
      if (update.has_value()) {
        if (is_ghost(v)) {
          auto ghost = v;
          buckets[ghosts.blocks[ghost - num_nodes()]].push_back(
              {.global_id = ghosts.global_ids[ghost - num_nodes()],
               .update = std::move(update.value())});
        } else {
          reached_blocks_by_border_node.clear();
          auto border_node = v;
          KASSERT(is_original_border_node(border_node));
          for (auto target : adjs[border_node]) {
            if (is_ghost(target)) {
              auto block = ghosts.blocks[target - num_nodes()];
              KASSERT(block < comm.size());
              if (reached_blocks_by_border_node.add(block)) {
                buckets[block].push_back({.global_id = border_node + get_first_global_node(),
                                          .update = std::move(update.value())});
              }
            }
          }
        }
      }
    }

    // return received ghost updates
    return with_flattened(buckets, comm.size()).call([&](auto... flattened) {
      return comm.alltoallv(std::move(flattened)...);
    });
  };

  template <typename QualifiedNodeID>
  std::vector<GlobalNodeID> ping_ghosts(sized_vector<QualifiedNodeID>& pings) const {
    using namespace kamping;

    std::vector<std::vector<GlobalNodeID>> buckets(comm.size());
    fast_set reached_blocks_by_border_node(comm.size());

    // scan updated border nodes and prepare update messages for each reached
    // blocks
    for (NodeID border_node : pings) {
      reached_blocks_by_border_node.clear();
      KASSERT(is_original_border_node(border_node));
      for (auto target : adjs[border_node]) {
        if (is_ghost(target)) {
          auto block = ghosts.blocks[target - num_nodes()];
          KASSERT(block < comm.size());
          if (reached_blocks_by_border_node.add(block)) {
            buckets[block].push_back(border_node + get_first_global_node());
          }
        }
      }
    }

    // return received ghost updates
    return with_flattened(buckets, comm.size()).call([&](auto... flattened) {
      return comm.alltoallv(std::move(flattened)...);
    });
  };

  const sized_vector<NodeID>& get_last_bulk_modified_nodes() { return bulk_buffer; }

  ghost_node_mapping take_ghosts_mapping() { return std::move(ghosts); }

  kamping::Communicator<> comm;  ////< kamping communicator (MPI)
 protected:
  GlobalNodeID first_global_node;  ////< first global node ID of distributed
                                   /// graph stored at this rank
  NodeID number_nodes;             ////< number of local nodes (ghosts excluded)
  NodeID visible_nodes;            ////< number of local visible nodes
  NodeID visible_ghosts;

  sized_vector<NodeID> border;  ////< border node IDs (at least one ghost as neighbor)
  fast_set border_set;

  // ghost nodes
  ghost_node_mapping ghosts;  ////< ghost mapper

  // access [0, num_nodes()+num_ghosts())
  std::vector<dynamic_neighborhood> adjs;  ////< dynamic neighborhoods
  std::vector<bool> hidden;                ////< maintain which nodes are visible
  sized_vector<NodeID> hidden_nodes;
  std::vector<NodeWeight> weights;  ////< maintain weights of nodes

  sized_vector<NodeID> bulk_buffer;
  fast_set bulk_buffer_set;
};

static_assert(std::is_move_constructible<dist_dynamic_graph>::value,
              "dist_dynamic_graph must have a move constructor!");

}  // namespace kadisredu
