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

using std::cout;
using std::endl;
using std::string;
using std::tuple;
using std::unordered_set;
using std::vector;

int main(int argc, char* argv[])
{
  // string directory = ARG_NEEDLE_TESTDATA_DIR "/length_5e4_samples_2e4/";

  // ARG arg = arg_utils::arg_from_ts_files(directory + "nodes.txt", directory + "edges.txt");
  // ARG arg = arg_utils::deserialize_arg_cpp("/gpfs3/well/palamara/users/ray826/arg_needle_lib_dev/test.argn");
  std::ifstream inputFile("/gpfs3/well/palamara/users/ray826/arg_needle_lib_dev/geno_test_loc_154401679.txt");

  std::vector<int> geno;
  std::string line;
  int geno_sum = 0;
  while (std::getline(inputFile, line)) {
      try {
          int number = std::stoi(line);
          geno.push_back(number);
          geno_sum += number;
      } catch (const std::invalid_argument& e) {
          std::cerr << "Invalid number in file: " << line << std::endl;
      } catch (const std::out_of_range& e) {
          std::cerr << "Number out of range in file: " << line << std::endl;
      }
  }

  inputFile.close();

  std::cout << "total carriers " << geno_sum << std::endl;

  ARG arg = arg_utils::deserialize_arg_cpp("/well/palamara/users/awo066/TDPBWT/experiments/huge_imputation_arg/results_v2/args/chr1.from143199864.to249222527.chunk1-2.argn", 154401678 - 144021214, 154401680 - 144021214);
  arg.populate_children_and_roots();
  DescendantList::set_threshold(64);
  // arg.check_basic(false);
  cout << arg.arg_nodes.size() << " nodes, " << arg.get_breakpoints().size() << " trees" << endl;
  cout << arg.arg_nodes.size() << " nodes, " << arg.num_edges() << " edges" << endl;

  arg_utils::map_genotype_to_ARG(arg, geno, 154401679 - arg.offset);

  // // auto bits = arg_utils::get_bitset(arg.arg_nodes.at(1296387).get(), arg.leaf_ids.size(), 339.85263061523438);
  // std::size_t k = 10;
  // Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic> matrix = Eigen::Matrix<int8_t, Eigen::Dynamic, Eigen::Dynamic>::Zero(k, arg.leaf_ids.size());

  // // Seed for random number generator
  // std::srand(std::time(0));

  // // Fill the matrix with random 0s and 1s
  // for (int i = 0; i < k; ++i) {
  //   for (int j = 0; j < arg.leaf_ids.size(); ++j) {
  //     if (std::rand() % 2000) matrix(i, j) = 0;
  //     else matrix(i, j) = 1;
  //   }
  // }

  // std::vector<double> positions(k);

  // for (int i = 0; i < k; ++i) {
  //   float randomValue = 0 + static_cast<double>(std::rand()) / (static_cast<double>(RAND_MAX / (5e4 - 0)));
  //   positions[i] = randomValue;
  // }

  // // Sort the vector
  // std::sort(positions.begin(), positions.end());

  // arg_utils::map_genotypes_to_ARG(arg, matrix, positions, 1);

  return 0;
}
