#ifndef PROPCA_ALSTRUCTURE_H_
#define PROPCA_ALSTRUCTURE_H_

#include <cassert>
#include <chrono>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <tuple>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SVD>

#include <Spectra/SymEigsSolver.h>

#include "arg.hpp"
#include "arg_edge.hpp"
#include "arg_node.hpp"
#include "arg_utils.hpp"
#include "descendant_list.hpp"
#include "random_utils.hpp"

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXdr;

struct ALSoptions {
  std::string GENOTYPE_FILE_PATH = "";
  std::string ROWSPACE_FILE_PATH = "";
  std::string INITIAL_FILE_PATH = "";
  std::string FREQ_FILE_PATH = "";
  std::string OUTPUT_PATH = "";
  int max_iterations = 1000;
  int num_of_evec = 5;
  bool debugmode = false;
  double convergence_limit = 0.00001;
  int seed = -1;
  int n_threads = 1;
  bool write_to_file = false;
};

class ALStructure {
public:
  using Scalar = double;
  ALSoptions opts;

  std::vector<ARG*> arg_list; // pointers to args passed in from python

  // ***********************************
  // *     ALStructure constructs      *
  // *  as defined in the SCOPE paper  *
  // ***********************************

  Eigen::VectorXd D; // (n,1) for LSE
  MatrixXdr V;        // (n,k) for truncated ALS
  MatrixXdr Fhat;     // (p,n) for truncated ALS
  MatrixXdr Phat;     // (p,k) for truncated ALS
  MatrixXdr Qhat;     // (k,n) for truncated ALS
  MatrixXdr Qhat_old; // (k,n) for truncated ALS
  MatrixXdr diff;     // (k,n) for truncated ALS
  double rmse;
  void calculate_D(void);
  void calculate_F(void);
  void solve_for_Qhat(void);
  void solve_for_Phat(void);
  void subspace_estimation(void);
  void initialize(std::default_random_engine& prng_eng);
  void truncated_alternating_least_squares(bool projection_mode = false);
  void fit_ALS();

  // ***********************************
  // *        PCA and spectral         *
  // *     clustering constructs       *
  // ***********************************

  Eigen::MatrixXd mean_allele_count; // (p,1) for PCA
  Eigen::VectorXd inv_var; // (p,1) inverse of variance for PCA
//   MatrixXdr centering_mat; // (p, p) for PCA
//   void fit_PCA(void);
//   void calculate_spectral_D(void);
//   void fit_spectral(void);
  

  // ***********************************
  // *        Shared variables         *
  // *        and IO functions         *
  // ***********************************

  int k, n; // number of eigenvectors in the algorithm
  unsigned p; // total number of mutations in (chunks of) ARG

  mutable unsigned int nops = 0u; // number of matrix-vector product operations
  long long int niter; // number of iterations in the solution
  double convergence_limit; // local copy from opts
  clock_t total_begin; // = clock();

  bool debug; // = false;

  std::vector<unsigned> p_list; // number of variants in each chunk
  std::vector<unsigned> n_list; // TODO: check this is a constant vector for chunks
  std::vector<int> p_split_positions; // Index of chunk splits
  std::string output_path;

  void write_matrix(MatrixXdr& mat, const std::string file_name);
  void write_matrix_maf(MatrixXdr& mat, const std::string file_name);
  void write_vector(Eigen::VectorXd& vec, const std::string file_name);

  // Genotype matrix-vector product 
  Eigen::MatrixXd parallel_post_multiplication(const Eigen::MatrixXd &in_mat) const;
  Eigen::MatrixXd parallel_pre_multiplication(const Eigen::MatrixXd &in_mat) const;
  unsigned int rows() const; // rows for Spectra
  unsigned int cols(); // cols for Spectra
  void perform_op(const double* x_in, double* y_out) const; // operation for Spectra

  ALStructure(std::vector<ARG*> input_arg_list, ALSoptions options = ALSoptions()) : arg_list(input_arg_list) {
    // Set default values
    arg_list = input_arg_list;
    nops = 1;
    opts = options;
  }
};

#endif // PROPCA_ALSTRUCTURE_H_