#include "arg_mult.hpp"

void ARGMatMult::load_arg(ARG& arg)
{
  if (arg.num_mutations() == 0) {
    THROW_LINE("No mutation found on the ARG.");
  }
  arg_utils::prepare_fast_multiplication(arg);

  std::unordered_map<int, std::map<double, std::set<int>>> node_to_split_to_mut_set_id;
  node_to_split_to_mut_set_id.reserve(arg.arg_nodes.size());
  // fill the split positions
  for (auto node_it = arg.fast_multiplication_data.topo_order.begin();
       node_it != arg.fast_multiplication_data.topo_order.end(); node_it++) {
    auto& splits = arg.fast_multiplication_data.node_id_to_split_points.at(*node_it);
    node_to_split_to_mut_set_id.emplace(*node_it, std::map<double, std::set<int>>());
    for (auto s : splits) {
      node_to_split_to_mut_set_id.at(*node_it).emplace(s, std::set<int>());
    }
  }

  // begin traversing the arg from root to leaves, filling in the topology of mutation sets
  int mut_set_id = 0;
  boost::timer::progress_display pbar(arg.num_mutations());

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
          std::set<int> current_split_result;
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
            node_to_split_to_mut_set_id.at(*node_it).at(current_split).merge(current_split_result);
          } else {
            std::set<int> current_split_muts;
            for (auto m : muts) {
              current_split_muts.emplace(arg.fast_multiplication_data.pos_to_mut_id.at(m->position));
              ++pbar;
            }
            mut_topo_desc_to_anc.emplace(mut_set_id, current_split_result);
            mut_set_id_to_muts.emplace(mut_set_id, current_split_muts);
            std::set<int> m_id;
            m_id.emplace(mut_set_id);
            node_to_split_to_mut_set_id.at(*node_it).at(current_split).merge(m_id);
            mut_set_id++;
          }
        }
      }
    }
  }

  // now fill the topology in reverse order, i.e. from leaf to root
  for (int i = 0; i != mut_set_id; i++)
    mut_topo_anc_to_desc.emplace(i, std::set<int>());

  for (auto& entry : mut_topo_desc_to_anc) {
    for (auto& anc : entry.second) {
      mut_topo_anc_to_desc.at(anc).emplace(entry.first);
    }
  }

  // also need a topological ordering of these mutation sets
  std::queue<int> topo_to_process; // sets with zero descendant
  std::unordered_map<int, int> desc_count;

  for (auto& mut_set_entry : mut_topo_anc_to_desc) {
    if (!mut_set_entry.second.empty()) {
      desc_count.emplace(mut_set_entry.first, mut_set_entry.second.size());
    } else {
      topo_to_process.push(mut_set_entry.first);
    }
  }

  mut_set_topo_order_leaf_to_root.reserve(mut_set_id);
  while (!topo_to_process.empty()) {
    int current_mut_set_id = topo_to_process.front();
    topo_to_process.pop();
    mut_set_topo_order_leaf_to_root.emplace_back(current_mut_set_id);
    auto& current_mut_set_anc = mut_topo_desc_to_anc.at(current_mut_set_id);
    for (auto& entry : current_mut_set_anc) {
      desc_count.at(entry) -= 1;
      if (desc_count.at(entry) == 0)
        topo_to_process.push(entry);
    }
  }
  desc_count.clear();

  // finally special treatment for the sample leaves, linking them with their immediate ancestral mutation sets
  for (auto& leaf_id : arg.leaf_ids) {
    std::set<int> muts;
    for (auto split : arg.fast_multiplication_data.node_id_to_split_points.at(leaf_id)) {
      muts.merge(node_to_split_to_mut_set_id.at(leaf_id).at(split));
    }
    indiv_to_mut_set_id.emplace(leaf_id, muts);
  }

  // save the count of mutations
  n_mut_indexed = arg.num_mutations();

  // save allele freq
  allele_frequencies = arg.fast_multiplication_data.allele_frequencies;

  // save leaf set size
  n_leaves = arg.leaf_ids.size();
}

Eigen::MatrixXd ARGMatMult::left_mult(
    const Eigen::MatrixXd& in_mat, bool standardize_mut, arg_real_t alpha, bool diploid)
{
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(in_mat.rows(), n_mut_indexed);

  // dimension check
  if (allele_frequencies.size() != n_mut_indexed) {
    throw std::runtime_error(THROW_LINE("Mismatching index data. Re-run ARGMatMult::load_arg(ARG& arg)."));
  }
  if (diploid) {
    if (n_leaves != in_mat.cols() * 2) {
      std::cout << "Leaf set size is " << n_leaves << " but vector size is 2 * " << in_mat.cols() << std::endl;
      throw std::runtime_error(THROW_LINE("Mismatching sample and vector sizes"));
    }
  } else {
    if (n_leaves != in_mat.cols()) {
      std::cout << "Leaf set size is " << n_leaves << " but vector size is " << in_mat.cols() << std::endl;
      throw std::runtime_error(THROW_LINE("Mismatching sample and vector sizes"));
    }
    if (standardize_mut) {
      throw std::runtime_error(THROW_LINE("Standardized haploid genotype is not yet implemented..."));
    }
  }

  // initialise the results starting from the leaves
  for (auto& entry : indiv_to_mut_set_id) {
    int indv_id = diploid ? entry.first / 2 : entry.first;
    auto target_col = in_mat.col(indv_id);
    for (auto anc_mut_set : entry.second) {
      for (int mut_id : mut_set_id_to_muts.at(anc_mut_set))
        result.col(mut_id) += target_col;
    }
  }

  // traverse upwards
  for (auto mut_set_id : mut_set_topo_order_leaf_to_root) {
    Eigen::VectorXd temp_result = Eigen::VectorXd::Zero(in_mat.rows());
    auto desc_mut_sets = mut_topo_anc_to_desc.at(mut_set_id);
    for (auto desc_mut_set_id : desc_mut_sets) {
      temp_result += result.col(*mut_set_id_to_muts.at(desc_mut_set_id).begin());
    }
    for (auto mut_id : mut_set_id_to_muts.at(mut_set_id)) {
      result.col(mut_id) += temp_result;
    }
  }

  // postprocessing and deal with normalisation
  if (standardize_mut) {
    Eigen::VectorXd mat_row_sum = in_mat.rowwise().sum();

    Eigen::ArrayXd means = allele_frequencies.array() / n_leaves;
    Eigen::VectorXd stds = (2 * means * (1 - means)).sqrt().pow(alpha);

    Eigen::MatrixXd offset = mat_row_sum.matrix() * (2*means).matrix().transpose();
    result -= offset;
    result *= stds.asDiagonal();
    // rowsum = input_mat.sum(axis=1)
    // af = geno.mean(axis=0, keepdims=True)/2
    // offset = rowsum.reshape(-1, 1) @ (2*af).reshape(1, -1)
    // std = np.sqrt(2 * af * (1-af))
    // return (input_mat @ geno - offset)/std
  }
  return result;
}

Eigen::MatrixXd ARGMatMult::right_mult(
    const Eigen::MatrixXd& in_mat, bool standardize_mut, arg_real_t alpha, bool diploid)
{

  // dimension check

  if (n_mut_indexed != in_mat.rows()) {
    std::cout << "Mutations are " << n_mut_indexed << " but vector size is " << in_mat.rows() << std::endl;
    throw std::runtime_error(THROW_LINE("Mismatching mutations and matrix row sizes"));
  }
  if (allele_frequencies.size() != n_mut_indexed) {
    throw std::runtime_error(THROW_LINE("Mismatching index data. Re-run ARGMatMult::load_arg(ARG& arg)."));
  }
  if (!diploid && standardize_mut) {
    throw std::runtime_error(
        THROW_LINE("Standardized haploid genotype is not yet implemented..."));
  }

  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(diploid ? n_leaves/2 : n_leaves, in_mat.cols());
  Eigen::MatrixXd input_copy = in_mat;
  Eigen::MatrixXd partial_results = Eigen::MatrixXd::Zero(mut_set_id_to_muts.size(), in_mat.cols());

  if (standardize_mut) {
    Eigen::VectorXd mat_row_sum = in_mat.rowwise().sum();
    Eigen::ArrayXd means = allele_frequencies.array() / n_leaves;
    Eigen::VectorXd stds = (2 * means * (1 - means)).sqrt().pow(alpha);
    input_copy = stds.asDiagonal() * input_copy;

  }

  // traverse downwards
  for (auto mut_set_id = mut_set_topo_order_leaf_to_root.rbegin(); mut_set_id != mut_set_topo_order_leaf_to_root.rend() ; mut_set_id++) {
    Eigen::VectorXd current_result = Eigen::VectorXd::Zero(in_mat.cols());
    auto anc_mut_sets = mut_topo_desc_to_anc.at(*mut_set_id);
    for (auto anc_mut_set_id : anc_mut_sets) {
      current_result += partial_results.row(anc_mut_set_id);
    }
    for (auto mut_id : mut_set_id_to_muts.at(*mut_set_id)) {
      current_result += input_copy.row(mut_id);
    }
    partial_results.row(*mut_set_id) += current_result;
  }

  // finish with the individuals
  for (auto& entry : indiv_to_mut_set_id) {
    for (auto mut_set_id : entry.second) {
      result.row(diploid ? entry.first / 2 : entry.first) += partial_results.row(mut_set_id);
    }
  }

  // postprocessing with normalisation
  if (standardize_mut) {
    Eigen::ArrayXd means = allele_frequencies.array() / n_leaves;
    Eigen::VectorXd stds = (2 * means * (1 - means)).sqrt().pow(alpha);
    Eigen::VectorXd offset = 2 * means.matrix().transpose() * stds.asDiagonal() * in_mat;
    result -= offset;
    // af = geno.mean(axis=0, keepdims=True)/2
    // std = np.sqrt(2 * af * (1-af))
    // offset = np.ones((n,1)) @ np.divide((2*af).reshape(1, -1), std.reshape(1,-1)) @ input_mat_rhs
    // return geno @ np.diag(1/std.flatten()) @ input_mat_rhs - offset
  }
  return result;
}