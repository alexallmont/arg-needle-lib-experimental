/*
  This file is part of the ARG-Needle genealogical inference and
  analysis software suite.
  Copyright (C) 2023 ARG-Needle Developers.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// Visit example ARGs in various ways

#include "arg.hpp"
#include "alstructure.hpp"
#include "arg_edge.hpp"
#include "arg_node.hpp"
#include "arg_utils.hpp"
#include "descendant_list.hpp"
#include "genotype_mapping.hpp"
#include "random_utils.hpp"
#include "serialize_arg.hpp"
#include "types.hpp"
#include "utils.hpp"


#include <cassert>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>
#include <fstream>
#include <boost/progress.hpp>

using std::cout;
using std::endl;
using std::string;
using std::tuple;
using std::unordered_set;
using std::vector;

int main(int argc, char* argv[])
{

  ARG arg = arg_utils::deserialize_arg_cpp("/well/palamara/users/awo066/TDPBWT/experiments/1K_genomes/EUR/results/chr1.from752721.to121475791/topology.argn");
  arg.populate_children_and_roots();
  arg_utils::generate_mutations(arg, 1e-9, 2);

  std::set<std::pair<int, arg_real_t>> keep_edge;

  for (auto& mut : arg.get_mutations()) {
    std::deque<ARGEdge*> edges_to_visit;
    edges_to_visit.push_front(mut->edge);

    while (!edges_to_visit.empty()) {
      auto current_edge = edges_to_visit.front();
      edges_to_visit.pop_front();
      keep_edge.emplace(current_edge->child->ID, current_edge->start);
      auto current_edge_children = current_edge->child->children_at(mut->position);
      if (!current_edge_children.empty()) {
        for (auto edge_ptr : current_edge_children) {
          edges_to_visit.push_back(edge_ptr);
        }
      }
    }
  }
  // DescendantList::set_threshold(64);
  // arg.check_basic(false);
  cout << arg.arg_nodes.size() << " nodes, " << arg.get_breakpoints().size() << " trees" << endl;
  cout << arg.leaf_ids.size() << " samples, " << arg.num_edges() << " edges" << endl;
  cout << keep_edge.size() << " edges are to be kept for " << arg.num_mutations() << " mutations" << endl;

  for (auto& node_it : arg.arg_nodes) {
    for (auto& edge_it : node_it.second->parents){

      if (keep_edge.find(std::make_pair(node_it.second->ID, edge_it.second->start)) == keep_edge.end()){
        node_it.second->remove_parent(edge_it.second->start);
      }
    }
  }
  arg.populate_children_and_roots();

  std::vector<int> nodes_to_remove;
  for (auto& node_it : arg.arg_nodes) {
    if ((node_it.second->children == nullptr || node_it.second->children->empty()) && !arg.is_leaf(node_it.second->ID)) {
      nodes_to_remove.push_back(node_it.second->ID);
    }
  }
  for (auto id : nodes_to_remove) {
    arg.arg_nodes.erase(id);
  }
  cout << arg.arg_nodes.size() << " nodes, " << arg.get_breakpoints().size() << " trees" << endl;
  cout << arg.leaf_ids.size() << " samples, " << arg.num_edges() << " edges" << endl;
  arg.populate_children_and_roots();
  


  arg_utils::prepare_fast_multiplication(arg);

  auto out = arg_utils::ARG_matrix_multiply_existing_mut_fast(arg, Eigen::MatrixXd::Ones(1, arg.leaf_ids.size()/2), false, 0, true);
  cout << "tot allele count " << out.sum() << endl;
  // std::vector<ARG*> arg_list;
  // arg_list.push_back(&arg);

  // auto als = ALStructure(arg_list);

  // als.fit_ALS();


  return 0;
}

