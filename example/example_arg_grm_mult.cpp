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
  // string directory = ARG_NEEDLE_TESTDATA_DIR "/length_5e3_samples_2e3/";

  // ARG arg = arg_utils::arg_from_ts_files(directory + "nodes.txt", directory + "edges.txt");
  // ARG arg = arg_utils::deserialize_arg_cpp("/gpfs3/well/palamara/users/ray826/arg_needle_lib_dev/test.argn");
  // std::ifstream inputFile("/gpfs3/well/palamara/users/ray826/arg_needle_lib_dev/geno_test_loc_154401679.txt");

  // std::vector<int> geno;
  // std::string line;
  // int geno_sum = 0;
  // while (std::getline(inputFile, line)) {
  //     try {
  //         int number = std::stoi(line);
  //         geno.push_back(number);
  //         geno_sum += number;
  //     } catch (const std::invalid_argument& e) {
  //         std::cerr << "Invalid number in file: " << line << std::endl;
  //     } catch (const std::out_of_range& e) {
  //         std::cerr << "Number out of range in file: " << line << std::endl;
  //     }
  // }
  // inputFile.close();

  // std::cout << "total carriers " << geno_sum << std::endl;

  ARG arg = arg_utils::deserialize_arg_cpp("/well/palamara/users/awo066/TDPBWT/experiments/huge_imputation_arg/results_v2/args/chr1.from143199864.to249222527.chunk1-2.argn");
  arg.populate_children_and_roots();
  // DescendantList::set_threshold(64);
  // arg.check_basic(false);
  cout << arg.arg_nodes.size() << " nodes, " << arg.get_breakpoints().size() << " trees" << endl;
  cout << arg.leaf_ids.size() << " samples, " << arg.num_edges() << " edges" << endl;

  arg_utils::prepare_fast_multiplication(arg);

  // int total_splits = 0;
  // for (auto& entry : arg.fast_multiplication_data.node_id_to_split_points) total_splits += entry.second.size();
  


  Eigen::MatrixXd input_mat = Eigen::MatrixXd::Ones(arg.leaf_ids.size(), 1);
  input_mat = Eigen::MatrixXd::Zero(arg.leaf_ids.size(), 1);
  input_mat(2,0) = 1;
  // map from (node_id, split_pos) -> partial results
  std::map< std::pair<int, arg_real_t>, Eigen::VectorXd> split_node_results; 

  // filling in place holder data
  for (auto entry : arg.fast_multiplication_data.node_id_to_split_points) {
    for (auto pos : entry.second) {
      int node_to_insert = entry.first;
      split_node_results[std::make_pair(node_to_insert, pos)] = Eigen::VectorXd::Zero(input_mat.cols());
    }
  }

  // initialize the leaves
  for (auto leaf_id : arg.leaf_ids) {
    split_node_results[std::make_pair(leaf_id, 0.)] = input_mat.row(leaf_id);
    // split_node_results[std::make_pair(leaf_id, arg.end)] = input_mat.row(leaf_id);
  }
  cout << split_node_results.size() << endl;


  boost::progress_display pbar(arg.arg_nodes.size());

  // first travel upwards from leaves to root
  for (auto it = arg.fast_multiplication_data.topo_order.rbegin(); it != arg.fast_multiplication_data.topo_order.rend(); it++) {
    auto node_id = *it;
    ++pbar;
    auto& node_parents = arg.arg_nodes.at(node_id)->parents;
    for (auto& parent_entry : node_parents) {
      auto parent_id = parent_entry.second.get()->parent->ID;
      auto edge_start = parent_entry.second.get()->start;
      auto edge_end = parent_entry.second.get()->end;

      // find the relevant split pos that spans the edge range
      auto &splits = arg.fast_multiplication_data.node_id_to_split_points.at(node_id);
      auto last_split = splits.lower_bound(edge_end);
      // need to backtrack one here since the ends are also included in node_id_to_split_points.
      // we only need to include the start positions in iteration below
      // last_split--;
      auto first_split = splits.upper_bound(edge_start);
      // again go back by one to find a starting split no later than supplied start_pos.
      first_split--;
      // cout << "node " << node_id << " with parent " << parent_id << " at edge " << edge_start << " to " << edge_end << endl;
      for (auto split_it = first_split; split_it != last_split; split_it++)
      {
        auto interval_left = *split_it;
        split_it++;
        auto interval_right = *split_it;
        split_it--;

        // cout << "   " << interval_left << " to " << interval_right << " "<< std::flush;

        // now update parent split results
        auto &parent_splits = arg.fast_multiplication_data.node_id_to_split_points.at(parent_id);
        auto first_parent_split = parent_splits.lower_bound(std::max(edge_start, interval_left));
        // no need to shift since children splits are always parent splits
        // assert(*first_parent_split == edge_start);
        auto last_parent_split = parent_splits.lower_bound(std::min(edge_end, interval_right));
        // if (last_parent_split == parent_splits.end()) last_parent_split--;
        assert(last_parent_split != parent_splits.end());
        for (auto it = first_parent_split; it != last_parent_split; it++) {
          split_node_results.at(std::make_pair(parent_id, *it)) += split_node_results.at(std::make_pair(node_id, interval_left));
          // cout << "(" << *it << "<-" << interval_left << ")" << std::flush;
          assert(interval_left <= *it);
        }

      }
      // cout << endl;

      // auto& parent_splits = arg.fast_multiplication_data.node_id_to_split_points.at(parent_id);
      // for (auto parent_split_it = parent_splits.lower_bound(edge_start); parent_split_it != parent_splits.lower_bound(edge_end); parent_split_it++) {
      //   auto parent_split_pos = *parent_split_it;

      //   // at this parent split point, the corresponding child descendants is recorded at the child split point exactly at or to the left of the split
      //   auto child_split_it = arg.fast_multiplication_data.node_id_to_split_points.at(node_id).upper_bound(parent_split_pos);
      //   child_split_it--;

      // }

    }
    
  }

  for (auto& entry : split_node_results) {
    // cout << entry.first.first << ", " << entry.first.second << " : " << entry.second << endl;
    assert(entry.second[0] <= arg.leaf_ids.size());
  }

  arg_real_t tot_vol = 0;
  Eigen::VectorXd tot_vec = Eigen::VectorXd::Zero(input_mat.cols());
  boost::progress_display pbar2(arg.arg_nodes.size());

  for (auto it = arg.fast_multiplication_data.topo_order.begin(); it != arg.fast_multiplication_data.topo_order.end(); it++) {
    auto node_id = *it;
    ++pbar2;
    auto& node_parents = arg.arg_nodes.at(node_id)->parents;

    auto split_it_start = arg.fast_multiplication_data.node_id_to_split_points.at(node_id).begin();
    auto split_it_end = arg.fast_multiplication_data.node_id_to_split_points.at(node_id).end();
    split_it_end--;

    for (auto split_it = split_it_start; split_it != split_it_end; split_it++) {
      auto split_pos = *split_it;
      split_it++;
      auto next_split_pos = *split_it;
      split_it--;
      auto& node_parents = arg.arg_nodes.at(node_id)->parents;
      Eigen::VectorXd current_cum_vector = Eigen::VectorXd::Zero(input_mat.cols());
      arg_real_t current_cum_vol = 0.;

      auto parent_edge_it_start = node_parents.upper_bound(split_pos);
      if (parent_edge_it_start != node_parents.begin()) parent_edge_it_start--;
      auto parent_edge_it_end = node_parents.lower_bound(next_split_pos);
      for (auto parent_edge_it = parent_edge_it_start; parent_edge_it!=parent_edge_it_end; parent_edge_it++) {
        auto edge_start = parent_edge_it->second->start;
        auto edge_end = parent_edge_it->second->end;
        if (edge_end <= split_pos) continue;
        auto parent_id = parent_edge_it->second->parent->ID;
        auto parent_split_lookup_start = std::max(edge_start, split_pos);
        auto parent_split_lookup_end = std::min(edge_end, next_split_pos);

        current_cum_vol += (parent_edge_it->second->parent->height - parent_edge_it->second->child->height)*(parent_split_lookup_end - parent_split_lookup_start);
        
        for (auto it = arg.fast_multiplication_data.node_id_to_split_points.at(parent_id).lower_bound(parent_split_lookup_start); it != arg.fast_multiplication_data.node_id_to_split_points.at(parent_id).lower_bound(parent_split_lookup_end); it++) {
          current_cum_vector += split_node_results.at(std::make_pair(parent_id, *it));
        }
      }

      split_node_results.at(std::make_pair(node_id, split_pos)) *= current_cum_vol;
      split_node_results.at(std::make_pair(node_id, split_pos)) += current_cum_vector;

      tot_vol += current_cum_vol;
      // cout << current_cum_vector << endl;
    }
  }


  // for (auto it = arg.fast_multiplication_data.topo_order.begin(); it != arg.fast_multiplication_data.topo_order.end(); it++) {
  //   auto node_id = *it;
  //   auto& node_parents = arg.arg_nodes.at(node_id)->parents;

  //   for (auto& split_pos : arg.fast_multiplication_data.node_id_to_split_points.at(node_id)) {
  //     if ((split_pos == arg.end) || (arg.arg_nodes.at(node_id)->parent_edge_at(split_pos) == nullptr))
  //     split_node_results.at(std::make_pair(node_id, split_pos)) = Eigen::VectorXd::Zero(input_mat.cols());
  //   }
  //   if (node_parents.empty()) {
  //     continue;
  //   }
  //   for (auto& parent_entry : node_parents) {
  //     auto parent_id = parent_entry.second.get()->parent->ID;
  //     auto edge_start = parent_entry.second.get()->start;
  //     auto edge_end = parent_entry.second.get()->end;

  //     // find the relevant split pos that spans the edge range
  //     auto &splits = arg.fast_multiplication_data.node_id_to_split_points.at(node_id);
  //     auto first_split = splits.upper_bound(edge_start);
  //     // again go back by one to find a starting split no later than supplied start_pos.
  //     first_split--;
  //     auto last_split = splits.lower_bound(edge_end);
  //     for (auto split_it = first_split; split_it != last_split; split_it++)
  //     {
  //       auto interval_left = *split_it;
  //       split_it++;
  //       auto interval_right = *split_it;
  //       split_it--;

  //       Eigen::VectorXd current_split_result = Eigen::VectorXd::Zero(input_mat.cols());

  //       // now update parent split results
  //       auto &parent_splits = arg.fast_multiplication_data.node_id_to_split_points.at(parent_id);
  //       auto first_parent_split = parent_splits.lower_bound(std::max(edge_start, interval_left));
  //       // no need to shift since children splits are always parent splits
  //       // assert(*first_parent_split == edge_start);
  //       auto last_parent_split = parent_splits.lower_bound(std::min(edge_end, interval_right));
  //       // if (last_parent_split == parent_splits.end()) last_parent_split--;
  //       assert(last_parent_split != parent_splits.end());
  //       for (auto it = first_parent_split; it != last_parent_split; it++) {

  //         current_split_result += split_node_results.at(std::make_pair(parent_id, *it));
  //         cout << "adding to node " << node_id << " at split "<< interval_left << " sum of " << split_node_results.at(std::make_pair(parent_id, *it))[0] << std::endl;
  //         assert(interval_left <= *it);
  //       }
  //       auto volume = (std::min(edge_end, interval_right) - std::max(edge_start, interval_left))*(parent_entry.second->parent->height - parent_entry.second->child->height);
  //       current_split_result += volume * split_node_results.at(std::make_pair(node_id, interval_left));

  //       cout << "local edge vol " << node_id << " at split "<< interval_left << " " << volume << " * " << split_node_results.at(std::make_pair(node_id, interval_left))[0] << std::endl;

  //       cout << split_node_results.at(std::make_pair(node_id, interval_left)) << endl;
  //       // assert(volume * split_node_results.at(std::make_pair(node_id, interval_left))[0] < 1e11);
  //       // cout << split_node_results.at(std::make_pair(node_id, interval_left)) << endl;
  //       // assert(current_split_result[0] < 1e11);
  //       std::cout << "assigning node " << node_id << " at split " << interval_left <<" "<< current_split_result << endl;
  //       cout << endl;
  //       split_node_results.at(std::make_pair(node_id, interval_left)) = current_split_result;

  //     }

  //   }
    
  // }
  cout << "done " << endl;
  std::cout << split_node_results.at(std::make_pair(0,0)) << std::endl;



  arg_real_t root_vol = 0;
  for (auto& entry : arg.roots) {
    root_vol += entry.second->node->height * (entry.second->end - entry.second->start);
  }
  cout << root_vol << endl;

  

  std::cout << tot_vol << std::endl;
  std::cout << arg_utils::total_volume(arg) << std::endl;

  std::cout << arg_utils::distance_matrix_v2(arg).at(0).at(1) << std::endl;
  // std::cout << arg_utils::distance_matrix(arg).at(0).at(0) << std::endl;

  return 0;
}
