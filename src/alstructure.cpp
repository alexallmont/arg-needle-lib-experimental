#include "alstructure.hpp"

using std::cout;
using std::endl;
using std::string;
using std::tuple;

// ################################
// ||                            ||
// ||  START OF CODE FROM SCOPE  ||
// ||                            ||
// ################################

struct timespec t0;

template <typename M> M load_tsv(const std::string& path);
template <typename M> M read_plink_freq_file(const std::string& path, int k);

double divide_by_two(double x) {
  return x / 2.0;
}

double truncate_with_epsilon(double x) {
  double epsilon = 0.0000000001;
  if (x <= 0.0)
    return epsilon;
  if (x >= 1.0)
    return 1.0 - epsilon;
  return x;
}

void project_onto_simplex(std::vector<double>& data) {
  std::vector<size_t> inds(data.size());

  std::iota(inds.begin(), inds.end(), 0);

  std::stable_sort(
      inds.begin(), inds.end(), [&data](size_t i, size_t j) { return data[i] > data[j]; });

  double tmpsum = 0;
  double tmax;
  bool bget = false;

  for (int i = 1; i < data.size(); i++) {
    tmpsum = tmpsum + data[inds[i - 1]];
    tmax = (tmpsum - 1.0) / i;
    if (tmax >= data[inds[i]]) {
      bget = true;
      break;
    }
  }

  if (!(bget)) {
    tmpsum = tmpsum + data[inds[data.size() - 1]];
    tmax = (tmpsum - 1) / data.size();
  }

  for (int i = 0; i < data.size(); i++) {
    data[i] = data[i] - tmax;
    if (data[i] < 0.0) {
      data[i] = 0.0;
    }
  }
}

void ALStructure::solve_for_Qhat() {
  MatrixXdr temp_kxp(k, p);
  temp_kxp = ((Phat.transpose() * Phat).inverse()) * Phat.transpose();
  // temp_kxn = temp_kxp * X, where X is the p X n genotype matrix
  MatrixXdr temp_kxn(k, n);
  temp_kxn = parallel_pre_multiplication(temp_kxp);
  temp_kxn = temp_kxn.unaryExpr(&divide_by_two);
  Qhat = (temp_kxn * V) * V.transpose();
}

void ALStructure::solve_for_Phat() {
  MatrixXdr temp_nxk(n, k);
  temp_nxk = V * (V.transpose() * ((Qhat.transpose() * (Qhat * Qhat.transpose()).inverse())));
  Phat = parallel_post_multiplication(temp_nxk);
  Phat = Phat.unaryExpr(&divide_by_two);
}

void ALStructure::initialize(std::default_random_engine& prng_eng) {
  if (opts.seed == -1) {
    opts.seed = static_cast<int>(time(NULL));
  }

  std::cout << "Initializing Phat using seed " << opts.seed << std::endl;
  prng_eng.seed(opts.seed);

  std::uniform_real_distribution<double> dis(0, 1);
  Phat = MatrixXdr::Zero(p, k).unaryExpr([&](float dummy) { return dis(prng_eng); });

  // std::cout << "Phat is populated with random values " << Phat << std::endl;

  if (opts.debugmode)
    write_matrix(Phat, std::string("Phat_0_") + std::to_string(opts.seed) + std::string(".txt"));
}

void ALStructure::truncated_alternating_least_squares(bool projection_mode) {
  if (projection_mode) {
    std::cout << "Solving for Q using provided frequencies" << std::endl;
    solve_for_Qhat();
    std::vector<double> col;
    col.resize(k);
    for (int c_iter = 0; c_iter < n; c_iter++) {
      // VectorXd::Map(&col[0], d) = Qhat.col(c_iter);
      for (int r_iter = 0; r_iter < k; r_iter++) {
        col[r_iter] = Qhat(r_iter, c_iter);
      }

      project_onto_simplex(col);

      for (int r_iter = 0; r_iter < k; r_iter++) {
        Qhat(r_iter, c_iter) = col[r_iter];
      }
    }
    return;
  }
  solve_for_Qhat();
  if (opts.debugmode)
    write_matrix(Qhat, std::string("Qhat_0.txt"));

  solve_for_Phat();
  if (opts.debugmode)
    write_matrix(Phat, std::string("Phat_1.txt"));

  for (niter = 1; niter < opts.max_iterations; niter++) {
    Qhat_old = Qhat;
    solve_for_Qhat();
    if (opts.debugmode)
      write_matrix(Qhat, std::string("Qhat_" + std::to_string(niter) + ".txt"));

    std::vector<double> col;
    col.resize(k);
    for (int c_iter = 0; c_iter < n; c_iter++) {
      // VectorXd::Map(&col[0], d) = Qhat.col(c_iter);
      for (int r_iter = 0; r_iter < k; r_iter++) {
        col[r_iter] = Qhat(r_iter, c_iter);
      }

      project_onto_simplex(col);

      for (int r_iter = 0; r_iter < k; r_iter++) {
        Qhat(r_iter, c_iter) = col[r_iter];
      }
    }
    if (opts.debugmode)
      write_matrix(Qhat, std::string("Qhat_" + std::to_string(niter) + "_w_constraints.txt"));

    solve_for_Phat();
    if (opts.debugmode)
      write_matrix(Phat, std::string("Phat_" + std::to_string(niter + 1) + ".txt"));

    Phat = Phat.unaryExpr(&truncate_with_epsilon);
    if (opts.debugmode)
      write_matrix(Phat, std::string("Phat_" + std::to_string(niter + 1) + "_w_constraints.txt"));

    diff = Qhat - Qhat_old;
    rmse = diff.norm() / sqrt(n * k);
    std::cout << "Iteration " << niter + 1 << "  -- RMSE " << std::setprecision(15) << rmse
              << std::endl;
    if ((rmse <= opts.convergence_limit) || std::isnan(rmse)) {
      std::cout << "Breaking after " << niter + 1 << " iterations" << std::endl;
      break;
    }
  }
}

void ALStructure::write_matrix(MatrixXdr& mat, const std::string file_name) {
  if (opts.write_to_file) {
    std::ofstream fp;
    fp.open((opts.OUTPUT_PATH + file_name).c_str());
    fp << std::setprecision(15) << mat << std::endl;
    fp.close();
  }
  else {
    std::cout << "Skipping file output" << std::endl;
  }
}

void ALStructure::write_matrix_maf(MatrixXdr& mat, const std::string file_name) {
  if (opts.write_to_file) {
    std::ofstream fp;
    fp.open((opts.OUTPUT_PATH + file_name).c_str());
    fp << std::setprecision(15) << 1 - mat.array() << std::endl;
    fp.close();
  }
  else {
    std::cout << "Skipping file output" << std::endl;
  }
}

void ALStructure::write_vector(Eigen::VectorXd& vec, const std::string file_name) {
  std::ofstream fp;
  fp.open((opts.OUTPUT_PATH + file_name).c_str());
  fp << std::setprecision(15) << vec << std::endl;
  fp.close();
}

void ALStructure::subspace_estimation() {
  std::cout << "Performing latent subspace estimation" << std::endl;

  std::cout << "Calculating D matrix" << std::endl;
  calculate_D();

  if (opts.debugmode)
    write_vector(D, "D.txt");

  // Calculate V
  std::cout << "Calculating V matrix" << std::endl;
  nops = 0;
  Spectra::SymEigsSolver<ALStructure> eigs(*this, k, k * 2 + 1);
  eigs.init();
  eigs.compute(Spectra::SortRule::LargestAlge, opts.max_iterations, opts.convergence_limit);

  if (eigs.info() == Spectra::CompInfo::Successful) {
    V = eigs.eigenvectors();
    write_matrix(V, "V.txt");
    if (opts.debugmode) {
      Eigen::VectorXd evals = eigs.eigenvalues().array() / (n - 1);
      write_vector(evals, "evals.txt");
    }
    std::cout << "Latent subspace esimation completed after " << nops << " iterations" << std::endl;
  }
  else {
    throw new std::runtime_error(std::string("Spectra eigendecomposition unsucessful") +
                                 ", status" + std::to_string(static_cast<int>(eigs.info())));
  }
}

// https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix
template <typename M> M load_tsv(const std::string& path) {
  std::ifstream indata;
  indata.open(path);
  std::string line;
  std::vector<double> values;
  int rows = 0;
  while (std::getline(indata, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    while (std::getline(lineStream, cell, '\t')) {
      values.push_back(std::stod(cell));
    }
    ++rows;
  }
  return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime,
                                        M::ColsAtCompileTime, Eigen::RowMajor>>(
      values.data(), rows, values.size() / rows);
}

template <typename M> M read_plink_freq_file(const std::string& path, int k) {
  /*
  Reads from plink.frq.strat file to initialize P matrix
  Returns rows/k x k matrix, where rows in number of rows
          in file
  */
  std::ifstream indata;
  indata.open(path);
  std::string line;
  std::vector<double> values;
  int rows = 1;
  std::getline(indata, line);
  while (std::getline(indata, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    std::vector<std::string> seglist;
    while (lineStream >> cell) {
      seglist.push_back(cell);
    }
    values.push_back(1 - std::stod(seglist[5]));
    ++rows;
  }
  return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime,
                                        M::ColsAtCompileTime, Eigen::RowMajor>>(
      values.data(), rows / k, k);
}

// ***********************************
// *           Fitting ALS           *
// ***********************************

void ALStructure::calculate_D() {
  D.resize(n);
  D.fill(0.);
  for (int i = 0; i < p_list.size(); i++) {
    Eigen::VectorXd D_chunk = Eigen::VectorXd::Zero(n);
    ARG* arg = arg_list[i];
    arg_utils::visit_mutations(
        *arg,
        [this, &D_chunk](DescendantList& desc_list, const Mutation* mutation) {
          // cout << mutation.height << endl;
          // we assert here as we don't know filtered #mutations when we start
          Eigen::VectorXd temp_pth_row_allel_count = Eigen::VectorXd::Zero(n);
          for (int v : desc_list.values()) {
            temp_pth_row_allel_count[v / 2] += 1;
          }
          // }
          D_chunk += 2 * temp_pth_row_allel_count -
                     temp_pth_row_allel_count.cwiseProduct(temp_pth_row_allel_count);
        },
        0, std::numeric_limits<arg_real_t>::infinity());
    { D += D_chunk; }
  }
}

void ALStructure::calculate_F() {
  // F is the allele frequency matrix as in ALStructure
  mean_allele_count.resize(p, 1);
  Eigen::MatrixXd in_mat(n, 1);
  Eigen::MatrixXd out_mat(p, 1);
  in_mat.fill(1.);
  out_mat = ALStructure::parallel_post_multiplication(in_mat);
  mean_allele_count = out_mat / 2 / n;
  inv_var.resize(p);
  Eigen::VectorXd temp(p);
  temp.fill(1.);
  temp -= mean_allele_count;           // (1-p)
  temp = 2 * mean_allele_count * temp; // 2*p*(1-p)
  inv_var = temp.cwiseInverse();
}

void ALStructure::fit_ALS() {

  total_begin = clock();
  clock_t io_begin = clock();

  int this_n_muts = 0;
  int this_n_haps = 0;
  n_list.clear();
  p_list.clear();

  for (auto& elem : arg_list) {
    this_n_muts = 0;
    p_list.push_back(elem->get_num_sites());

    this_n_haps = elem->leaf_ids.size();
    n_list.push_back(this_n_haps / 2);
  }
  p = std::accumulate(p_list.begin(), p_list.end(), 0);
  n = this_n_haps / 2;
  k = opts.num_of_evec;
  p_split_positions.resize(p_list.size());
  std::partial_sum(p_list.begin(), p_list.end(), p_split_positions.begin());

  auto start = std::chrono::system_clock::now();

  clock_gettime(CLOCK_REALTIME, &t0);

  clock_t io_end = clock();

  std::cout << "Running on Dataset of " << n << " individuals containing " << p
            << " SNPs" << std::endl;

  clock_t it_begin = clock();

  if (std::string(opts.ROWSPACE_FILE_PATH) != "") {
    std::cout << "Using provided V" << std::endl;

    // Read eigenvectors of the n x n matrix: G = (1/m) * (X^T X - D)
    V = load_tsv<MatrixXdr>(opts.ROWSPACE_FILE_PATH);
    if (k != V.cols()) {
      k = V.cols();
      std::cout << "Mismatch between column number of provided V and provided k!" << std::endl;
      std::cout << "Changing k to number of columns in V" << std::endl;
    }

    if (V.rows() != n) {
      std::cout << "Dimensions of genotype matrix and rowspace matrix do not agree!" << std::endl;
      exit(-1);
    }
  }
  else {
    subspace_estimation();
  }

  // Create pseudo random number generator (PRNG) engine
  std::default_random_engine prng_eng{};

  if (std::string(opts.INITIAL_FILE_PATH) != "") {
    std::cout << "Using initial Phat provided" << std::endl;
    Phat = load_tsv<MatrixXdr>(opts.INITIAL_FILE_PATH);
  }
  else if (std::string(opts.FREQ_FILE_PATH) != "") {
    Phat = read_plink_freq_file<MatrixXdr>(opts.FREQ_FILE_PATH, k);
    opts.max_iterations = 1;
  }
  else {
    initialize(prng_eng);
  }

  if (std::string(opts.FREQ_FILE_PATH) != "") {
    truncated_alternating_least_squares(true);
  }
  else {
    truncated_alternating_least_squares();

    // Try restarting with new seed!
    if (std::isnan(rmse) && (std::string(opts.INITIAL_FILE_PATH) == "")) {
      // opts.given_seed = false;
      for (int xx = 0; xx < 5; xx++) {
        initialize(prng_eng);
        truncated_alternating_least_squares();
        if (!std::isnan(rmse)) {
          break;
        }
      }
    }
  }

  clock_t it_end = clock();

  write_matrix_maf(Phat, "Phat.txt");
  write_matrix(Qhat, "Qhat.txt");

  clock_t total_end = clock();
  double io_time = static_cast<double>(io_end - io_begin) / CLOCKS_PER_SEC;
  // double avg_it_time = static_cast<double>(it_end - it_begin) / (MAX_ITER * 1.0 *
  // CLOCKS_PER_SEC);
  double total_time = static_cast<double>(total_end - total_begin) / CLOCKS_PER_SEC;
  std::cout << "Completed!" << std::endl;
  std::cout << "IO Time:  " << io_time << std::endl;
  std::cout << "Total runtime:   " << total_time << std::endl;

  std::chrono::duration<double> wctduration = std::chrono::system_clock::now() - start;
  std::cout << "Wall clock time = " << wctduration.count() << std::endl;
}

// ################################
// ||                            ||
// ||   END OF CODE FROM SCOPE   ||
// ||                            ||
// ################################



// ***********************************
// *      Matrix-vector product      *
// ***********************************

unsigned int ALStructure::cols() {
  return n;
}

unsigned int ALStructure::rows() const {
  return n;
}

void ALStructure::perform_op(const double* x_in, double* y_out) const {
    // Performs ((Xv)^T X)^T - Dv
    MatrixXdr x = Eigen::Map<const Eigen::VectorXd>(x_in, n);
    Eigen::Map<Eigen::VectorXd> y(y_out, n);

    MatrixXdr temp_px1(p, 1);

    temp_px1 = parallel_post_multiplication(x);
    temp_px1.transposeInPlace(); // now 1xp
    MatrixXdr temp_nx1(n, 1);

    temp_nx1 = parallel_pre_multiplication(temp_px1);
    temp_nx1.transposeInPlace(); // now nx1
    y.noalias() = temp_nx1 - D.cwiseProduct(x);
    nops++;
}

Eigen::MatrixXd ALStructure::parallel_post_multiplication(const Eigen::MatrixXd& in_mat) const {
  assert(in_mat.rows() == n);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(in_mat.cols(), p);
  {
    for (int i = 0; i < p_list.size(); i++) {
      int p_block_begin = (i == 0) ? 0 : p_split_positions[i - 1];
      MatrixXdr p_chunk(p_list[i], in_mat.cols());
      ARG* arg = arg_list[i];
      p_chunk = arg_utils::ARG_matrix_multiply_existing_mut_fast_mt(
          *arg, in_mat.transpose(), false, 0, true, opts.n_threads);
      {
        result.block(0, p_block_begin, in_mat.cols(), p_list[i]) = p_chunk;
        // std::cout << "thread " << omp_get_thread_num() << " is invoked in post mult" <<
        // std::endl;
      }
    }
  }
  result.transposeInPlace();
  return result;
}

Eigen::MatrixXd ALStructure::parallel_pre_multiplication(const Eigen::MatrixXd& in_mat) const {
  assert(in_mat.cols() == p);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, in_mat.rows());
  {
    for (int i = 0; i < p_list.size(); i++) {
      int p_block_begin = (i == 0) ? 0 : p_split_positions[i - 1];
      ARG* arg = arg_list[i];
      Eigen::MatrixXd in_mat_chunk(in_mat.rows(), p_list[i]);
      in_mat_chunk = in_mat.block(0, p_block_begin, in_mat.rows(), p_list[i]);
      in_mat_chunk.transposeInPlace();
      result += arg_utils::ARG_matrix_multiply_samples_faster_mt(
          *arg, in_mat_chunk, false, 0, true, opts.n_threads);
    }
  }
  result.transposeInPlace();
  return result;
}