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

#include "alstructure.hpp"
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

#include <boost/functional/hash.hpp>
#include <boost/timer/progress_display.hpp>
#include <cassert>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>

using std::cout;
using std::endl;
using std::pair;
using std::set;
using std::stack;
using std::string;
using std::tuple;
using std::unordered_set;
using std::vector;

int main(int argc, char* argv[])
{

  // ARG arg =
  // arg_utils::deserialize_arg_cpp("/well/palamara/users/awo066/TDPBWT/experiments/1K_genomes/EUR/results/chr1.from752721.to121475791/topology.argn");
  ARG arg = arg_utils::deserialize_arg_cpp("/well/palamara/users/ray826/arg_needle_lib_dev/test_giant.argn");
  arg.populate_children_and_roots();
  arg_utils::generate_mutations(arg, 1e-8, 22);
  arg_utils::prepare_fast_multiplication(arg);
  boost::timer::progress_display pbar(arg.num_mutations());

  std::unordered_map<int, std::set<int>> mut_topo;
  std::unordered_map<std::pair<int, double>, std::set<int>, boost::hash<std::pair<int, double>>> node_split_pos_to_muts;
  std::map<int, std::map<double, set<int>>> node_to_split_to_mut_set_id;
  std::map<int, set<int>> mut_set_id_to_muts;

  for (auto node_it = arg.fast_multiplication_data.topo_order.begin();
       node_it != arg.fast_multiplication_data.topo_order.end(); node_it++) {
    auto& splits = arg.fast_multiplication_data.node_id_to_split_points.at(*node_it);
    node_to_split_to_mut_set_id.emplace(*node_it, std::map<double, set<int>>());
    for (auto s : splits) {
      node_to_split_to_mut_set_id.at(*node_it).emplace(s, set<int>());
    }
  }


  int mut_set_id=0;
  for (auto node_it = arg.fast_multiplication_data.topo_order.begin();
       node_it != arg.fast_multiplication_data.topo_order.end(); node_it++) {
    auto& splits = arg.fast_multiplication_data.node_id_to_split_points.at(*node_it);
    auto& parents = arg.arg_nodes.at(*node_it)->parents;

    if (!parents.empty()) {
      for (auto& parent_edge_entry : parents) {
        auto p_start = parent_edge_entry.first;
        auto& p_edge = parent_edge_entry.second;
        auto split_start = splits.upper_bound(p_start);
        split_start--;
        auto split_end = splits.lower_bound(p_edge->end);
        for (auto split_it = split_start; split_it != split_end; split_it++) {
          auto current_split = *split_it;
          split_it++;
          auto next_split = *split_it;
          split_it--;
          auto sum_start = std::max(p_edge->start, current_split);
          auto sum_end = std::min(p_edge->end, next_split);
          // we sum up results for each split by counting contributions from each parent edge
          set<int> current_split_result;
          auto parent_res_st = node_to_split_to_mut_set_id.at(p_edge->parent->ID).lower_bound(sum_start);
          assert(sum_start == parent_res_st->first);
          auto parent_res_ed = node_to_split_to_mut_set_id.at(p_edge->parent->ID).lower_bound(sum_end);
          assert(sum_end == parent_res_ed->first);
          for (auto parent_res = parent_res_st; parent_res != parent_res_ed; parent_res++) {
            auto e = parent_res->second;
            current_split_result.merge(e);
          }
          auto muts = p_edge->mutations_in_range(sum_start, sum_end);
          if (muts.empty()) {
            // if (node_to_split_to_muts.at(*node_it).find(sum_start) == node_to_split_to_muts.at(*node_it).end()) {
            //   node_to_split_to_muts.at(*node_it).emplace(sum_start, current_split_result);
            // }
            // else {
            node_to_split_to_mut_set_id.at(*node_it).at(current_split).merge(current_split_result);
            // }
          } else {
            set<int> current_split_muts;
            for (auto m : muts) {
              current_split_muts.emplace(arg.fast_multiplication_data.pos_to_mut_id.at(m->position));
              // mut_topo.emplace(arg.fast_multiplication_data.pos_to_mut_id.at(m->position), current_split_result);
              ++pbar;
              // cout << "mut at pos " << m->position << " id " <<
              // arg.fast_multiplication_data.pos_to_mut_id.at(m->position) << " is desc to mut ids "; for (auto k :
              // current_split_result) {
              //   cout << k << " ";
              // }
              // cout << endl;
            }
            // if (node_to_split_to_muts.at(*node_it).find(sum_start) == node_to_split_to_muts.at(*node_it).end()) {
            //   node_to_split_to_muts.at(*node_it).emplace(sum_start, current_split_muts);
            // }
            // else {
            mut_topo.emplace(mut_set_id, current_split_result);
            mut_set_id_to_muts.emplace(mut_set_id, current_split_muts);
            set<int> m_id;
            m_id.emplace(mut_set_id);
            node_to_split_to_mut_set_id.at(*node_it).at(current_split).merge(m_id);
            mut_set_id++;

            // }
          }
        }
      }
    }
  }
  std::map<int, set<int>> indiv_to_mut;
  for (auto& arg_node_entry : arg.arg_nodes) {
    if (arg.is_leaf(arg_node_entry.first)) {
      set<int> muts;
      for (auto split : arg.fast_multiplication_data.node_id_to_split_points.at(arg_node_entry.first)) {
        muts.merge(node_to_split_to_mut_set_id.at(arg_node_entry.first).at(split));
      }
      indiv_to_mut.emplace(arg_node_entry.first, muts);
    }
  }

  Eigen::MatrixXi mut_mat = arg_utils::get_mutations_matrix(arg).cast<int>(); // p x n
  // cout << mut_mat.rows() << " " << mut_mat.cols() << endl;

  // for (auto& entry : mut_topo) {
  //   auto desc_mut = entry.first;
  //   // auto desc_row = mut_mat.row(desc_mut);
  //   cout << "checking mut set id " << desc_mut << " including mut id {";
  //   for (auto m : mut_set_id_to_muts.at(desc_mut)) {
  //     cout << m << " ";
  //     assert(mut_mat.row(m) == mut_mat.row(*mut_set_id_to_muts.at(desc_mut).begin()));
  //   }
  //   cout << "}\ndescendant to\n";
  //   if (!entry.second.empty()){
  //     for (auto& anc : entry.second) {
  //       cout << "  mut set id " << anc << " including mut id {" << std::flush;
  //       for (auto m : mut_set_id_to_muts.at(anc)) {
  //         cout << m << " " << std::flush;
  //         assert(mut_mat.row(m) == mut_mat.row(*mut_set_id_to_muts.at(anc).begin()));
  //       }
  //       cout << "}\n";
  //       // cout << " " << anc;
  //     }

  //   }
  // }

  cout << mut_set_id << endl;
  cout << arg.num_mutations() << endl;
  cout << arg.fast_multiplication_data.pos_to_mut_id.size() << endl;

  Eigen::VectorXi found_mut(arg.num_mutations());
  found_mut.fill(0);
  for (auto& entry : mut_set_id_to_muts) {
    for (auto& mid : entry.second) {
      assert(found_mut[mid] == 0);
      found_mut[mid] = 1;
    }
  }
  assert(found_mut.sum() == arg.num_mutations());

  std::map<int, set<int>> mut_topo_bottom_up;
  for (int i = 0; i != mut_set_id; i++)
    mut_topo_bottom_up.emplace(i, set<int>());

  for (auto& entry : mut_topo) {
    for (auto& anc : entry.second) {
      mut_topo_bottom_up.at(anc).emplace(entry.first);
    }
  }
  for (auto& entry : indiv_to_mut) {
    for (auto& mut : entry.second) {
      mut_topo_bottom_up.at(mut).emplace(-1 * (1+entry.first));
    }
  }


  cout << "verify genotypes" << endl;
  boost::timer::progress_display pgbar(mut_set_id);

  for (auto& entry : mut_topo_bottom_up) {
    ++pgbar;
    auto anc_mut = entry.first;
    // cout << "checking mut set id " << anc_mut << " including mut id { " << std::flush;
    Eigen::RowVectorXi actual_geno = mut_mat.row(*mut_set_id_to_muts.at(anc_mut).begin()).cast<int>();
    Eigen::RowVectorXi pred_geno = Eigen::VectorXi::Zero(arg.leaf_ids.size());
    // for (auto m : mut_set_id_to_muts.at(anc_mut)) {
    //   cout << m << " ";
    // }
    // cout << "}\n";
    if (!entry.second.empty()){
      for (auto& des : entry.second) {
        // cout << "  mut set id " << des << " including mut id {" << std::flush;
        if (des >= 0) {
          // for (auto m : mut_set_id_to_muts.at(des)) {
          //   cout << m << " " << std::flush;
          // }
          // cout << "}\n";
          pred_geno += mut_mat.row(*mut_set_id_to_muts.at(des).begin()).cast<int>();
        }
        else {
          // cout << "indv_" << (-1*des - 1) << "}" << endl;
          pred_geno[-1*des - 1] += 1;
        }
      }

    }
    if(actual_geno != pred_geno) {
      cout << actual_geno << endl;
      cout << pred_geno << endl;
      break;
    }
  }
  cout << "verify indv genotypes" << endl;
  boost::timer::progress_display ppbar(arg.leaf_ids.size());

  for (auto& leaf_id : arg.leaf_ids) {
    ++ppbar;
    // cout << "checking mut set id " << anc_mut << " including mut id { " << std::flush;
    Eigen::RowVectorXi actual_geno = mut_mat.col(leaf_id).cast<int>();
    Eigen::RowVectorXi pred_geno = Eigen::VectorXi::Zero(arg.num_mutations());

    std::deque<int> mut_set_id_to_process;

    for (auto entry : indiv_to_mut.at(leaf_id)) {
      mut_set_id_to_process.emplace_back(entry);
    }

    while (!mut_set_id_to_process.empty()) {
      int mut_set_id = mut_set_id_to_process.front();
      mut_set_id_to_process.pop_front();
      for (int mut_id : mut_set_id_to_muts.at(mut_set_id)) {
        pred_geno[mut_id] += 1;
      }
      for (int more_mut_set : mut_topo.at(mut_set_id)) mut_set_id_to_process.emplace_back(more_mut_set);
    }
    if(actual_geno != pred_geno) {
      cout << actual_geno << endl;
      cout << pred_geno << endl;
      break;
    }
  }

  // for (auto& entry : mut_topo_bottom_up) {
  //   Eigen::RowVectorXi geno = Eigen::VectorXi::Zero(arg.leaf_ids.size());
  //   Eigen::RowVectorXi ref_geno = mut_mat.row(entry.first);
  //   cout << "check anc id " << entry.first;
  //   for (auto& desc : entry.second) {
  //     cout << " " << desc;
  //     if (desc < 0)
  //       geno[-1 * desc - 1] = 1; // singleton
  //     else {
  //       // cout << "\n+" << mut_mat.row(desc) << endl;
  //       geno += mut_mat.row(desc);
  //       for (int i = 0; i < geno.size(); i++) {
  //         if (geno[i] > 0)
  //           geno[i] = 1;
  //       }
  //     }
  //   }
  //   cout << endl;
  //   // cout << geno << endl;
  //   // cout << ref_geno << endl;
  //   assert(geno == ref_geno);
  // }
  // cout << mut_mat.row(1001) << endl;
  // cout << mut_mat.row(875) << endl;
  // cout << mut_mat.row(905) << endl;
  // cout << mut_mat.row(1049) << endl;
  // cout << mut_mat.row(17) << endl;
  // cout << mut_mat.row(34) << endl;

  // for (auto node_it = arg.fast_multiplication_data.topo_order.begin(); node_it !=
  // arg.fast_multiplication_data.topo_order.end(); node_it++) {
  //   auto& splits = arg.fast_multiplication_data.node_id_to_split_points.at(*node_it);
  //   set<double> split_with_mut;
  //   cout << *node_it << endl;
  //   // cout << arg.arg_nodes.at(*node_it)->parents.empty() << endl;
  //   if (!arg.arg_nodes.at(*node_it)->parents.empty()){
  //     for (auto& parent_edge_entry : arg.arg_nodes.at(*node_it)->parents) {
  //       if (parent_edge_entry.second->mutations != nullptr) {
  //         for (auto* mut : *parent_edge_entry.second->mutations){
  //           auto mut_pos = mut->position;
  //           auto split_for_mut = splits.upper_bound(mut_pos);
  //           split_for_mut--;
  //           split_with_mut.emplace(*split_for_mut);
  //         }
  //       }
  //     }
  //   }

  //   for (auto& split_pos : splits) {
  //     if (split_pos == *splits.rbegin()) continue;
  //     auto* parent_edge = arg.arg_nodes.at(*node_it)->parent_edge_at(split_pos);
  //     if (split_with_mut.find(split_pos) != split_with_mut.end()) {
  //       double next_split_pos = *splits.upper_bound(split_pos);
  //       for (auto mut : *parent_edge->mutations) {
  //         if (mut->position >= split_pos && mut->position < next_split_pos)
  //         {
  //           cout << *node_it << " has mut at " << mut->position << " between " << split_pos << " and " <<
  //           next_split_pos << endl; if (node_split_pos_to_muts.find(std::make_pair(*node_it, split_pos)) ==
  //           node_split_pos_to_muts.end()) {
  //             node_split_pos_to_muts.emplace(std::make_pair(*node_it, split_pos), std::set<int>());
  //           }
  //           node_split_pos_to_muts.at(std::make_pair(*node_it,
  //           split_pos)).emplace(arg.fast_multiplication_data.pos_to_mut_id.at(mut->position));

  //           //also sum up the results and write out
  //           auto parent_split_st =
  //           arg.fast_multiplication_data.node_id_to_split_points.at(parent_edge->parent->ID).lower_bound(split_pos);
  //           auto parent_split_ed =
  //           arg.fast_multiplication_data.node_id_to_split_points.at(parent_edge->parent->ID).lower_bound(next_split_pos);
  //           set<int> current_res;
  //           for (auto it=parent_split_st; it != parent_split_ed; it++) {
  //             auto p = node_split_pos_to_muts.at(std::make_pair(parent_edge->parent->ID, *it));
  //             current_res.merge(p);
  //           }
  //           mut_topo.emplace(arg.fast_multiplication_data.pos_to_mut_id.at(mut->position), current_res);
  //         }
  //       }
  //     }
  //     else {
  //       // the case of no mut, so copy the mut info downwards
  //       if (parent_edge == nullptr) {
  //         // this is when we start with the root so it has no parents
  //         if (node_split_pos_to_muts.find(std::make_pair(*node_it, split_pos)) == node_split_pos_to_muts.end()) {
  //             node_split_pos_to_muts.emplace(std::make_pair(*node_it, split_pos), std::set<int>());
  //           }
  //         else THROW_LINE("??? root node somehow has been filled with some parent data");
  //       }
  //       else {
  //         // this is when the node simply has an empty parent edge
  //         double next_split_pos = *splits.upper_bound(split_pos);
  //         auto parent_split_st =
  //         arg.fast_multiplication_data.node_id_to_split_points.at(parent_edge->parent->ID).lower_bound(split_pos);
  //         auto parent_split_ed =
  //         arg.fast_multiplication_data.node_id_to_split_points.at(parent_edge->parent->ID).lower_bound(next_split_pos);
  //         set<int> current_res;
  //         for (auto it=parent_split_st; it != parent_split_ed; it++) {
  //           cout << "here" << endl;
  //           cout << "parents splits in " << split_pos << " -> " << next_split_pos << endl;
  //           cout << "parent node " << parent_edge->parent->ID << "  split at " << *it << endl;
  //           cout << "max split pos " <<
  //           *arg.fast_multiplication_data.node_id_to_split_points.at(parent_edge->parent->ID).rbegin() << endl; cout
  //           << "seek end " << *parent_split_ed << endl; auto p =
  //           node_split_pos_to_muts.at(std::make_pair(parent_edge->parent->ID, *it)); current_res.merge(p); cout <<
  //           "pass" << endl;
  //         }
  //         node_split_pos_to_muts.emplace(std::make_pair(*node_it, split_pos), current_res);
  //       }
  //     }
  //   }
  // }

  return 0;
}
