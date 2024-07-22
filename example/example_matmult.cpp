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
#include "serialize_arg.hpp"
#include "descendant_list.hpp"
#include "random_utils.hpp"
#include "types.hpp"
#include "utils.hpp"

#include <cassert>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>

using std::cout;
using std::endl;
using std::string;
using std::tuple;
using std::unordered_set;
using std::vector;

int main(int argc, char* argv[]) {
  // string directory = ARG_NEEDLE_TESTDATA_DIR "/length_5e4_samples_2e4";
  // if (argc > 1) {
  //   directory = (string) argv[1];
  // }
  // if (directory[directory.size() - 1] != '/') {
  //   directory += "/";
  // }

  // if (argc > 2) {
  //   string s = argv[2];
  //   int threshold = atoi(s.c_str());
  //   cout << "Setting threshold to " << threshold << endl << endl;
  //   DescendantList::set_threshold(threshold);
  // }

  // // bool verbose = false;
  // // if (argc > 3 && (string) argv[3] != "0") {
  // //   verbose = true;
  // // }

  // ARG arg = arg_utils::arg_from_ts_files(directory + "nodes.txt", directory + "edges.txt");
  ARG arg = arg_utils::deserialize_arg_cpp("/gpfs3/well/palamara/users/ray826/arg_needle_lib_dev/test.argn");
  arg.populate_children_and_roots();
  arg.check_basic(false);
  cout << arg.arg_nodes.size() << " nodes, " << arg.get_breakpoints().size() << " trees" << endl;
  cout << arg.arg_nodes.size() << " nodes, " << arg.num_edges() << " edges" << endl;

  arg_utils::generate_mutations(arg, 1e-8, 18);

  cout << "generated " << arg.num_mutations() << " mutations" << endl;
  arg.populate_mutations_on_edges();

  arg_utils::prepare_fast_multiplication(arg);

  Eigen::MatrixXd random_input = Eigen::MatrixXd::Random(arg.num_mutations(), 500);
  Eigen::MatrixXd random_input2 = Eigen::MatrixXd::Random(500, arg.leaf_ids.size()/2);

  auto result = arg_utils::ARG_matrix_multiply_samples_faster(arg, random_input, true, -1);
  auto result2 = arg_utils::ARG_matrix_multiply_existing_mut_fast(arg, random_input2, true, -1);

  // auto geno = arg_utils::get_mutations_matrix(arg);

  return 0;
}
