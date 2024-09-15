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
#include <chrono>


using std::cout;
using std::endl;
using std::string;
using std::tuple;
using std::unordered_set;
using std::vector;

int main(int argc, char* argv[]) {
  // ARG arg = arg_utils::deserialize_arg_cpp("/gpfs3/well/palamara/users/ray826/fastmult_arg/single-chunk-test/args/neGBR_n1e+05_l30e6.argn");
  auto t1 = std::chrono::high_resolution_clock::now();
  ARG arg = arg_utils::deserialize_arg_cpp("/gpfs3/well/palamara/users/ray826/fastmult_arg/geno-grm-generated/args/chr4.from52684820.to190906015.chunk7-8.split.0.0.to.12261607.0.argn");

  auto t2 = std::chrono::high_resolution_clock::now();
  auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  std::cout << "read in took " << time_span.count() << " seconds." << endl;

  t1 = std::chrono::high_resolution_clock::now();
  arg.populate_children_and_roots();

  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  std::cout << "pop children took " << time_span.count() << " seconds." << endl;

  t1 = std::chrono::high_resolution_clock::now();
  arg_utils::generate_mutations(arg, 1e-10, 18);

  cout << "generated " << arg.num_mutations() << " mutations" << endl;

  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  std::cout << "gen muts took " << time_span.count() << " seconds." << endl;

  t1 = std::chrono::high_resolution_clock::now();
  arg_utils::prepare_fast_multiplication(arg);

  t2 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
  std::cout << "prep mult took " << time_span.count() << " seconds." << endl;
  // arg.populate_mutations_on_edges();

  // Eigen::MatrixXd rand_pheno = Eigen::MatrixXd::Random(arg.leaf_ids.size()/2, 20);

  // // auto result = arg_utils::association_mutation_fast(arg, rand_pheno);
  // auto result = arg_utils::weighted_mut_squared_norm(arg, Eigen::MatrixXd::Random(arg.leaf_ids.size()/2, 7), true);


  return 0;
}
