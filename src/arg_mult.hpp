#include <Eigen/Core>
#include <Eigen/Dense>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/functional/hash.hpp>
#include <boost/timer/progress_display.hpp>
#include <cassert>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>

#include "arg.hpp"
#include "arg_edge.hpp"
#include "arg_node.hpp"
#include "arg_utils.hpp"
#include "descendant_list.hpp"
#include "random_utils.hpp"
#include "H5Cpp.h"

class ARGMatMult
{
public:

  // convenient constructor to load the arg directly
  ARGMatMult(ARG& arg);

  // default constructor
  ARGMatMult() = default;

  // mapping from a descendant set of mutations to all of its ancestral sets of mutations
  std::vector<std::vector<int>> mut_topo_desc_to_anc;
  // mapping from an ancestral set of mutations to all of its descendant sets of mutations
  std::vector<std::vector<int>> mut_topo_anc_to_desc;
  // mapping from mutation set id to its member mutations
  std::vector<std::vector<int>> mut_set_id_to_muts;
  // mapping from leaf_id to the mutation set ids it contains
  std::vector<std::vector<int>> indiv_to_mut_set_id;
  // mapping from mutation set id to leaf ids that it is ancestral to
  std::vector<std::vector<int>> mut_set_id_to_indiv;
  // number of mutations indexed, for now it's the same as arg.num_mutations()
  int n_mut_indexed;
  // topological ordering of mutation sets
  std::vector<int> mut_set_topo_order_leaf_to_root;
  // allele freq
  Eigen::VectorXd allele_frequencies;
  // number of individuals
  int n_leaves;


  // build the skeleton topology from an arg with mutation
  void load_arg(ARG& arg);

  // multiply (n x p) genotype matrix from the left hand side, output is (k x p)
  Eigen::MatrixXd left_mult(const Eigen::MatrixXd& in_mat, bool standardize_mut=true, arg_real_t alpha=-1, bool diploid=true);
  // multiply (n x p) genotype matrix from the right hand side, output is (n x k)
  Eigen::MatrixXd right_mult(const Eigen::MatrixXd& in_mat, bool standardize_mut=true, arg_real_t alpha=-1, bool diploid=true);


  void save_hdf5(const std::string &filename, bool compress = false) ;
  void load_hdf5(const std::string &filename) ;

};