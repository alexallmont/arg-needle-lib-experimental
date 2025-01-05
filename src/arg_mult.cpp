#include "arg_mult.hpp"

ARGMatMult::ARGMatMult(ARG& arg) {
  ARGMatMult::load_arg(arg);
}

void ARGMatMult::load_arg(ARG& arg)
{
  if (arg.num_mutations() == 0) {
    THROW_LINE("No mutation found on the ARG.");
  }
  arg_utils::prepare_fast_multiplication(arg);

  boost::unordered_flat_map<int, boost::container::flat_map<double, std::vector<int>>> node_to_split_to_mut_set_id;
  node_to_split_to_mut_set_id.reserve(arg.arg_nodes.size());
  // fill the split positions
  for (auto node_it = arg.fast_multiplication_data.topo_order.begin();
       node_it != arg.fast_multiplication_data.topo_order.end(); node_it++) {
    auto& splits = arg.fast_multiplication_data.node_id_to_split_points.at(*node_it);
    node_to_split_to_mut_set_id.emplace(*node_it, boost::container::flat_map<double, std::vector<int>>());
    for (auto s : splits) {
      node_to_split_to_mut_set_id.at(*node_it).emplace(s, std::vector<int>());
    }
  }

  // begin traversing the arg from root to leaves, filling in the topology of mutation sets
  int mut_set_id = 0;
  mut_topo_anc_to_desc.reserve(arg.num_mutations());
  mut_topo_desc_to_anc.reserve(arg.num_mutations());
  mut_set_id_to_muts.reserve(arg.num_mutations());

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
          std::vector<int> current_split_result;
          auto parent_res_st = node_to_split_to_mut_set_id.at(p_edge->parent->ID).lower_bound(sum_start);
          assert(sum_start == parent_res_st->first);
          auto parent_res_ed = node_to_split_to_mut_set_id.at(p_edge->parent->ID).lower_bound(sum_end);
          assert(sum_end == parent_res_ed->first);
          for (auto parent_res = parent_res_st; parent_res != parent_res_ed; parent_res++) {
            auto e = parent_res->second;
            current_split_result.reserve(current_split_result.size() + std::distance(e.begin(), e.end()));
            current_split_result.insert(std::end(current_split_result), std::begin(e), std::end(e));
          }
          auto muts = p_edge->mutations_in_range(sum_start, sum_end);
          if (muts.empty()) {
            auto& current_node_result = node_to_split_to_mut_set_id.at(*node_it).at(current_split);
            current_node_result.reserve(current_node_result.size() + std::distance(current_split_result.begin(), current_split_result.end()));
            current_node_result.insert(std::end(current_node_result), std::begin(current_split_result), std::end(current_split_result));
            // node_to_split_to_mut_set_id.at(*node_it).at(current_split).merge(current_split_result);
          } else {
            std::vector<int> current_split_muts;
            for (auto m : muts) {
              current_split_muts.emplace_back(arg.fast_multiplication_data.pos_to_mut_id.at(m->position));
              ++pbar;
            }
            mut_topo_desc_to_anc.emplace_back(current_split_result);
            mut_set_id_to_muts.emplace_back(current_split_muts);
            std::vector<int> m_id;
            m_id.emplace_back(mut_set_id);
            auto& current_node_result = node_to_split_to_mut_set_id.at(*node_it).at(current_split);
            current_node_result.reserve(current_node_result.size() + std::distance(m_id.begin(), m_id.end()));
            current_node_result.insert(std::end(current_node_result), std::begin(m_id), std::end(m_id));
            // node_to_split_to_mut_set_id.at(*node_it).at(current_split).merge(m_id);
            mut_set_id++;
          }
        }
      }
    }
  }

  // now fill the topology in reverse order, i.e. from leaf to root
  for (int i = 0; i != mut_set_id; i++)
    mut_topo_anc_to_desc.emplace_back(std::vector<int>());

  for (int i=0; i != mut_set_id; i++) {
    for (auto& anc : mut_topo_desc_to_anc[i]) {
      mut_topo_anc_to_desc[anc].emplace_back(i);
    }
  }

  // also need a topological ordering of these mutation sets
  std::queue<int> topo_to_process; // sets with zero descendant
  std::vector<int> desc_count;
  desc_count.reserve(mut_set_id);

  for (int i=0; i != mut_set_id; i++) {
    int dc = mut_topo_anc_to_desc[i].size();
    desc_count.emplace_back(dc);
    if (dc == 0) {
      topo_to_process.push(i);
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
  for (int leaf_id=0; leaf_id != arg.leaf_ids.size(); leaf_id++) {
    std::vector<int> muts;
    for (auto split : arg.fast_multiplication_data.node_id_to_split_points.at(leaf_id)) {
      auto& mut_sets = node_to_split_to_mut_set_id.at(leaf_id).at(split);
      muts.reserve(muts.size() + std::distance(mut_sets.begin(), mut_sets.end()));
      muts.insert(std::end(muts), std::begin(mut_sets), std::end(mut_sets));
    }
    indiv_to_mut_set_id.emplace_back(muts);
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
  Eigen::MatrixXd partial_result = Eigen::MatrixXd::Zero(in_mat.rows(), mut_set_id_to_muts.size());

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

  auto t1 = std::chrono::high_resolution_clock::now();
  // initialise the results starting from the leaves
  // for (const auto& entry : indiv_to_mut_set_id) {
  //   int indv_id = diploid ? entry.first / 2 : entry.first;
  //   const auto& target_col = in_mat.col(indv_id);
  //   for (int anc_mut_set : entry.second) {
  //     // for (int mut_id : mut_set_id_to_muts.at(anc_mut_set))
  //     partial_result.col(anc_mut_set) += target_col;
  //   }
  // }
  for (int leaf_id=0; leaf_id != indiv_to_mut_set_id.size(); leaf_id++) {
    int indv_id = diploid ? leaf_id / 2 : leaf_id;
    const auto& target_col = in_mat.col(indv_id);
    for (int anc_mut_set : indiv_to_mut_set_id[leaf_id]) {
      // for (int mut_id : mut_set_id_to_muts.at(anc_mut_set))
      partial_result.col(anc_mut_set) += target_col;
    }
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  std::cout << "init build time " << duration.count() / 1000.0 << " s" << std::endl;

  t1 = std::chrono::high_resolution_clock::now();
  // traverse upwards
  for (auto mut_set_id : mut_set_topo_order_leaf_to_root) {
    Eigen::VectorXd temp_result = Eigen::VectorXd::Zero(in_mat.rows());
    auto desc_mut_sets = mut_topo_anc_to_desc.at(mut_set_id);
    for (auto desc_mut_set_id : desc_mut_sets) {
      // temp_result += result.col(*mut_set_id_to_muts.at(desc_mut_set_id).begin());
      temp_result += partial_result.col(desc_mut_set_id);
    }
    // for (auto mut_id : mut_set_id_to_muts.at(mut_set_id)) {
    //   result.col(mut_id) += temp_result;
    // }
    partial_result.col(mut_set_id) += temp_result;
  }
  t2 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  std::cout << "traverse time " << duration.count() / 1000.0 << " s" << std::endl;

  t1 = std::chrono::high_resolution_clock::now();
  for (int i=0; i!= mut_set_id_to_muts.size(); i++) {
    auto res = partial_result.col(i);
    for (int m_id : mut_set_id_to_muts[i]) {
      result.col(m_id) = res;
    }
  }
  t2 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  std::cout << "assign time " << duration.count() / 1000.0 << " s" << std::endl;

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

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> result = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(diploid ? n_leaves/2 : n_leaves, in_mat.cols());
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> input_copy = in_mat;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> partial_results = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(mut_set_id_to_muts.size(), in_mat.cols());

  if (standardize_mut) {
    Eigen::VectorXd mat_row_sum = in_mat.rowwise().sum();
    Eigen::ArrayXd means = allele_frequencies.array() / n_leaves;
    Eigen::VectorXd stds = (2 * means * (1 - means)).sqrt().pow(alpha);
    input_copy = stds.asDiagonal() * input_copy;

  }

  auto t1 = std::chrono::high_resolution_clock::now();
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
  auto t2 = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  std::cout << "traverse time " << duration.count() / 1000.0 << " s" << std::endl;

  t1 = std::chrono::high_resolution_clock::now();
  // finish with the individuals
  // for (auto& entry : indiv_to_mut_set_id) {
  //   for (auto mut_set_id : entry.second) {
  //     result.row(diploid ? entry.first / 2 : entry.first) += partial_results.row(mut_set_id);
  //   }
  // }
  for (int i=0; i!=indiv_to_mut_set_id.size(); i++) {
    for (auto mut_set_id : indiv_to_mut_set_id[i]) {
      result.row(diploid ? i / 2 : i) += partial_results.row(mut_set_id);
    }
  }
  t2 = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  std::cout << "assign time " << duration.count() / 1000.0 << " s" << std::endl;

  // postprocessing with normalisation
  if (standardize_mut) {
    Eigen::ArrayXd means = allele_frequencies.array() / n_leaves;
    Eigen::VectorXd stds = (2 * means * (1 - means)).sqrt().pow(alpha);
    Eigen::VectorXd offset = 2 * means.matrix().transpose() * stds.asDiagonal() * in_mat;
    result.rowwise() -= offset.transpose();
    // af = geno.mean(axis=0, keepdims=True)/2
    // std = np.sqrt(2 * af * (1-af))
    // offset = np.ones((n,1)) @ np.divide((2*af).reshape(1, -1), std.reshape(1,-1)) @ input_mat_rhs
    // return geno @ np.diag(1/std.flatten()) @ input_mat_rhs - offset
  }
  return result;
}