/*
  This file is part of the ARG-Needle genealogical inference and
  analysis software suite.
  Copyright (C) 2023-2025 ARG-Needle Developers.

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

#include "serialize_arg.hpp"

#include "deserialization_params.hpp"

#include "H5Cpp.h"

#include <iostream>
#include <filesystem>
#include <utility>
#include <vector>

using std::vector;

namespace
{

bool check_attribute(const H5::H5File& h5_file, const std::string& expected_attr)
{
  if (!h5_file.attrExists(expected_attr)) {
    std::cerr << "Expected file " << h5_file.getFileName() << " to include attribute `" << expected_attr << "`"
              << std::endl;
    return false;
  }
  return true;
}

bool check_dataset(const H5::H5File& h5_file, const std::string& expected_dset)
{
  try {
    auto dset = h5_file.openDataSet(expected_dset);
    return true;
  } catch (const H5::Exception&) {
    std::cerr << "Expected file " << h5_file.getFileName() << " to include dataset `" << expected_dset << "`"
              << std::endl;
    return false;
  }
}

bool check_group(const H5::H5File& h5_file, const std::string& expected_group)
{
  try {
    auto group = h5_file.openGroup(expected_group);
    return true;
  } catch (const H5::Exception&) {
    std::cerr << "Expected file " << h5_file.getFileName() << " to include group `" << expected_group << "`"
              << std::endl;
    return false;
  }
}

bool read_bool_attribute(const H5::H5File& file, const std::string& attrName)
{
  bool value = false;

  try {
    // Open the attribute from the file
    const H5::Attribute attribute = file.openAttribute(attrName);

    // Read the attribute, assuming it was stored as H5T_NATIVE_HBOOL
    uint8_t buffer{};
    attribute.read(attribute.getDataType(), &buffer);

    // Convert the integer buffer to bool
    value = (buffer != 0u);
  } catch (const H5::Exception& e) {
    // Handle exceptions: attribute not found, wrong type, etc.
    std::cerr << "Error reading attribute `" << attrName << "`: " << e.getDetailMsg() << std::endl;
  }

  return value;
}

template <class T> struct dependent_false : std::false_type {
};

template <class T> inline constexpr bool dependent_false_v = dependent_false<T>::value;

template <typename T> H5::PredType getPredType()
{
  if constexpr (std::is_same_v<T, double>) {
    return H5::PredType::NATIVE_DOUBLE;
  } else if constexpr (std::is_same_v<T, int>) {
    return H5::PredType::NATIVE_INT;
  } else if constexpr (std::is_same_v<T, uint8_t>) {
    return H5::PredType::NATIVE_UINT8;
  } else {
    static_assert(dependent_false_v<T>, "Unsupported type");
  }
}

// Templated function to read a dataset into a std::vector<T>
template<typename T>
std::vector<T> read_dataset_to_vector_1d(const H5::H5File& h5_file, const std::string& dset_name, hssize_t start = 0, hssize_t stop = -1) {
  std::vector<T> data;

  try {
    H5::DataSet dataset = h5_file.openDataSet(dset_name);
    H5::DataSpace dataspace = dataset.getSpace();

    if (dataspace.getSimpleExtentNdims() != 1) {
      throw std::runtime_error("Dataset must be 1-dimensional");
    }

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims);
    const hsize_t num_elements = dims[0];

    // Adjust stop value if necessary
    if (stop == -1 || stop > num_elements) {
      stop = num_elements;
    }

    if (start >= stop) {
      throw std::runtime_error("Invalid range: start must be less than stop");
    }

    const hsize_t block_size = stop - start;
    data.resize(block_size);

    // Define hyperslab in the dataset
    hsize_t offset[1] = {static_cast<hsize_t>(start)};
    hsize_t count[1] = {block_size};
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

    // Define the memory dataspace
    H5::DataSpace memspace(1, count);
    dataset.read(data.data(), getPredType<T>(), memspace, dataspace);
  } catch (const H5::Exception& e) {
    throw std::runtime_error("Failed to read dataset `" + dset_name + "`: " + e.getDetailMsg());
  }

  return data;
}

template <typename T>
std::vector<std::array<T, 2>> read_dataset_to_vector_2d(
    const H5::H5File& h5_file, const std::string& dset_name, hssize_t start = 0, hssize_t stop = -1)
{
  std::vector<std::array<T, 2>> data;

  try {
    H5::DataSet dataset = h5_file.openDataSet(dset_name);
    H5::DataSpace dataspace = dataset.getSpace();

    if (dataspace.getSimpleExtentNdims() != 2) {
      throw std::runtime_error("Dataset must be 2-dimensional");
    }

    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims);
    if (dims[1] != 2) {
      throw std::runtime_error("Second dimension of the dataset must be 2");
    }

    hsize_t num_elements = dims[0];

    // Adjust stop value if necessary
    if (stop == -1 || stop > num_elements) {
      stop = num_elements;
    }

    if (start >= stop) {
      throw std::runtime_error("Invalid range: start must be less than stop");
    }

    hsize_t block_size = stop - start;
    data.resize(block_size);

    // Define hyperslab in the dataset
    hsize_t offset[2] = {static_cast<hsize_t>(start), 0}; // Start from 'start', covering the entire 2nd dimension
    hsize_t count[2] = {block_size, 2}; // Number of elements to read in the 1st dim and include all of the 2nd dim
    dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

    // Define the memory dataspace
    H5::DataSpace memspace(2, count);
    dataset.read(data.data(), getPredType<T>(), memspace, dataspace);
  } catch (const H5::Exception& e) {
    throw std::runtime_error("Failed to read dataset `" + dset_name + "`: " + e.getDetailMsg());
  }

  return data;
}

int read_int_attribute(const H5::H5File& file, const std::string& attrName)
{
  int value{};

  try {
    // Open the attribute from the file
    H5::Attribute attribute = file.openAttribute(attrName);

    // Read the attribute, assuming it was stored as H5T_NATIVE_INT
    attribute.read(H5::PredType::NATIVE_INT, &value);
  } catch (const H5::Exception& e) {
    // Handle exceptions: attribute not found, wrong type, etc.
    std::cerr << "Error reading attribute `" << attrName << "`: " << e.getDetailMsg() << std::endl;
    exit(0);
  }

  return value;
}

bool validate_serialized_arg_v1(const H5::H5File& h5_file)
{
  // Expected attributes and datasets
  std::vector<std::string> expected_attrs = {"num_nodes", "num_edges", "num_mutations", "offset", "chromosome",
      "sequence_length", "datetime_created", "arg_file_version"};
  std::vector<std::string> expected_dsets = {"flags", "times", "edge_ranges", "edge_ids"};

  bool is_valid = true;

  for (const auto& attr : expected_attrs) {
    is_valid = check_attribute(h5_file, attr);
  }

  for (const auto& dset : expected_dsets) {
    is_valid = check_dataset(h5_file, dset);
  }

  return is_valid;
}

bool validate_serialized_arg_v2(const H5::H5File& h5_file)
{
  // Expected attributes and datasets
  std::vector<std::string> expected_attrs = {"num_nodes", "num_edges", "node_bounds", "num_mutations", "mutations",
      "offset", "chromosome", "start", "end", "threaded_samples", "datetime_created", "arg_file_version"};

  std::vector<std::string> expected_dsets = {"flags", "times", "edge_ranges", "edge_ids"};
  std::vector<std::string> optional_dsets = {"node_bounds"};

  std::vector<std::string> expected_groups = {};
  std::vector<std::string> optional_groups = {"mutations"};

  bool is_valid = true;

  for (const auto& attr : expected_attrs) {
    is_valid = check_attribute(h5_file, attr);
  }

  // The existence of optional datasets is marked by bool attributes of the same name
  for (const auto& dset_name : optional_dsets) {
    if (read_bool_attribute(h5_file, dset_name)) {
      expected_dsets.emplace_back(dset_name);
    }
  }

  for (const auto& dset : expected_dsets) {
    is_valid = check_dataset(h5_file, dset);
  }

  // The existence of optional groups is also marked by bool attributes of the same name
  for (const auto& group_name : optional_groups) {
    if (read_bool_attribute(h5_file, group_name)) {
      expected_groups.emplace_back(group_name);
    }
  }

  for (const auto& group : expected_groups) {
    is_valid = check_group(h5_file, group);
  }

  return is_valid;
}

ARG deserialize_arg_v1(const H5::H5File& h5_file, const int reserved_samples)
{

  // Read attributes
  const int offset = read_int_attribute(h5_file, "offset");
  const int chromosome = read_int_attribute(h5_file, "chromosome");
  const int sequence_length = read_int_attribute(h5_file, "sequence_length");
  const int num_nodes = read_int_attribute(h5_file, "num_nodes");
  const int num_edges = read_int_attribute(h5_file, "num_edges");

  // Read datasets
  std::vector<uint8_t> raw_flags = read_dataset_to_vector_1d<uint8_t>(h5_file, "flags");
  std::vector<double> raw_times = read_dataset_to_vector_1d<double>(h5_file, "times");
  std::vector<int> raw_edge_ids = read_dataset_to_vector_1d<int>(h5_file, "edge_ids");
  std::vector<double> raw_edge_ranges = read_dataset_to_vector_1d<double>(h5_file, "edge_ranges");

  assert(raw_flags.size() == num_nodes);
  assert(raw_times.size() == num_nodes);
  assert(raw_edge_ids.size() == 2 * num_edges);
  assert(raw_edge_ranges.size() == 2 * num_edges);

  // Translate from raw data to required data structures
  std::deque<bool> is_sample;
  for (const auto flag : raw_flags) {
    is_sample.push_back(flag != 0);
  }

  std::vector<std::pair<int, int>> edge_ids;
  for (std::size_t i = 0; i < raw_edge_ids.size(); i += 2) {
    edge_ids.emplace_back(raw_edge_ids[i], raw_edge_ids[i + 1]);
  }

  std::vector<std::pair<double, double>> edge_ranges;
  for (std::size_t i = 0; i < raw_edge_ranges.size(); i += 2) {
    edge_ranges.emplace_back(raw_edge_ranges[i], raw_edge_ranges[i + 1]);
  }

  // Construct the ARG object
  ARG arg(0, sequence_length, raw_times, is_sample, edge_ids, edge_ranges, reserved_samples);
  arg.set_offset(offset);
  arg.set_chromosome(chromosome);

  return arg;
}

ARG deserialize_arg_v2(const H5::H5File& h5_file, const int chunk_size, const int reserved_samples)
{
  DeserializationParams dp;
  dp.start = read_int_attribute(h5_file, "start");
  dp.end = read_int_attribute(h5_file, "end");
  dp.num_nodes = read_int_attribute(h5_file, "num_nodes");
  dp.offset = read_int_attribute(h5_file, "offset");
  dp.chromosome = read_int_attribute(h5_file, "chromosome");
  dp.threaded_samples = read_int_attribute(h5_file, "threaded_samples");
  dp.reserved_samples = reserved_samples;

  ARG arg(dp);

  // Process {chunk_size} nodes at a time, adding each chunk to the ARG as we go
  {
    const auto num_nodes = static_cast<hssize_t>(dp.num_nodes);
    hssize_t num_nodes_written = 0;

    while (num_nodes_written < num_nodes) {
      const hssize_t range_lo = num_nodes_written;
      const hssize_t range_hi = std::min(num_nodes_written + chunk_size, num_nodes);

      const auto node_heights = read_dataset_to_vector_1d<double>(h5_file, "times", range_lo, range_hi);
      const auto is_sample = read_dataset_to_vector_1d<uint8_t>(h5_file, "flags", range_lo, range_hi);

      if (read_bool_attribute(h5_file, "node_bounds")) {
        const auto node_bounds_data = read_dataset_to_vector_2d<double>(h5_file, "node_bounds", range_lo, range_hi);
        arg.deserialize_add_nodes(node_heights, is_sample, node_bounds_data);
      } else {
        arg.deserialize_add_nodes(node_heights, is_sample);
      }

      num_nodes_written += static_cast<hssize_t>(node_heights.size());
    }
  }

  // Process {chunk_size} edges at a time, adding each chunk to the ARG as we go
  {
    const auto num_edges = static_cast<hssize_t>(read_int_attribute(h5_file, "num_edges"));
    hssize_t num_edges_written = 0;

    while (num_edges_written < num_edges) {

      const hssize_t range_lo = num_edges_written;
      const hssize_t range_hi = std::min(num_edges_written + chunk_size, num_edges);

      const auto edge_id_data = read_dataset_to_vector_2d<int>(h5_file, "edge_ids", range_lo, range_hi);
      const auto edge_range_data = read_dataset_to_vector_2d<double>(h5_file, "edge_ranges", range_lo, range_hi);

      arg.deserialize_add_edges(edge_id_data, edge_range_data);

      num_edges_written += static_cast<hssize_t>(edge_id_data.size());
    }
  }

  // Process {chunk_size} mutations at a time, adding each chunk to the ARG as we go
  if (read_bool_attribute(h5_file, "mutations")) {

    const auto num_mutations = static_cast<hssize_t>(read_int_attribute(h5_file, "num_mutations"));
    hssize_t num_mutations_written = 0;

    while (num_mutations_written < num_mutations) {

      const hssize_t range_lo = num_mutations_written;
      const hssize_t range_hi = std::min(num_mutations_written + chunk_size, num_mutations);

      const auto mut_pos = read_dataset_to_vector_1d<double>(h5_file, "mutations/positions", range_lo, range_hi);
      const auto mut_hts = read_dataset_to_vector_1d<double>(h5_file, "mutations/heights", range_lo, range_hi);
      const auto mut_sid = read_dataset_to_vector_1d<int>(h5_file, "mutations/site_ids", range_lo, range_hi);
      const auto mut_eid = read_dataset_to_vector_2d<int>(h5_file, "mutations/edge_ids", range_lo, range_hi);

      arg.deserialize_add_mutations(mut_pos, mut_hts, mut_sid, mut_eid);

      num_mutations_written += static_cast<hssize_t>(mut_pos.size());
    }
  }

  return arg;
}

} // namespace

bool arg_utils::validate_serialized_arg(const std::string& file_name)
{
  // Check if file exists
  if (!std::filesystem::exists(file_name)) {
    std::cout << "File: " << file_name << " is not a valid file" << std::endl;
    return false;
  }

  if (!H5Fis_hdf5(file_name.c_str()))
  {
    std::cout << "File: " << file_name << " is not a valid HDF5 file" << std::endl;
    return false;
  }

  try {
    H5::H5File h5_file(file_name, H5F_ACC_RDONLY);

    // Check for the 'arg_file_version' attribute
    if (!check_attribute(h5_file, "arg_file_version")) {
      std::cout << "File: " << file_name
                << " is not a valid arg file because it does not contain `arg_file_version` attribute" << std::endl;
      return false;
    }

    const int arg_file_version = read_int_attribute(h5_file, "arg_file_version");

    // Validate file version
    if (arg_file_version == 1) {
      return validate_serialized_arg_v1(h5_file);
    }
    if (arg_file_version == 2) {
      return validate_serialized_arg_v2(h5_file);
    }

    std::cout << "Arg file version (" << arg_file_version << ") is not supported; valid versions are 1, 2."
              << std::endl;
    return false;

  } catch (const H5::Exception& e) {
    std::cerr << "HDF5 error on file: " << file_name << std::endl;
    std::cerr << e.getDetailMsg() << std::endl;
    return false;
  }
}

ARG arg_utils::deserialize_arg(const std::string& file_name, const int chunk_size, const int reserved_samples)
{
  if (!validate_serialized_arg(file_name)) {
    throw std::runtime_error("Invalid ARG file: " + file_name);
  }

  try {
    H5::H5File h5_file(file_name, H5F_ACC_RDONLY);

    const int arg_file_version = read_int_attribute(h5_file, "arg_file_version");

    if (arg_file_version == 1) {
      return deserialize_arg_v1(h5_file, reserved_samples);
    }
    if (arg_file_version == 2) {
      return deserialize_arg_v2(h5_file, chunk_size, reserved_samples);
    }

    throw std::logic_error(
        "Reached an unsupported arg_file_version after validation: " + std::to_string(arg_file_version));

  } catch (const H5::Exception& e) {
    std::cerr << "HDF5 error on file: " << file_name << std::endl;
    std::cerr << e.getDetailMsg() << std::endl;
    throw std::runtime_error("Unable to deserialize arg file: " + file_name);
  }


}

ARG arg_utils::deserialize_arg_cpp(const std::string& file_path, const arg_real_t trim_start, const arg_real_t trim_end, const arg_real_t truncation_height) {

  if (!std::filesystem::is_regular_file(file_path)) {
    throw std::logic_error("ARG file_path is not a path to a file");
  }

  // Try block to detect exceptions raised by any of the calls inside it
  try {
    // Turn off the auto-printing when failure occurs so that we can
    // handle the errors appropriately
    H5::Exception::dontPrint();

    // Open an existing file and dataset.
    H5::H5File arg_archive_file(file_path, H5F_ACC_RDONLY);

    int attr_arg_file_version{}, attr_chromosome{}, attr_num_edges{}, attr_num_mutations{},
        attr_offset{}, attr_num_nodes{};
    arg_real_t attr_sequence_length{};
    arg_archive_file.openAttribute("arg_file_version")
        .read(H5::PredType::NATIVE_INT, &attr_arg_file_version);

    if (attr_arg_file_version == 2) {
      // std::cerr << "Version 2 is not yet fully supported. \n";
    }

    arg_archive_file.openAttribute("num_nodes").read(H5::PredType::NATIVE_INT, &attr_num_nodes);
    arg_archive_file.openAttribute("chromosome").read(H5::PredType::NATIVE_INT, &attr_chromosome);
    arg_archive_file.openAttribute("num_edges").read(H5::PredType::NATIVE_INT, &attr_num_edges);
    arg_archive_file.openAttribute("num_mutations")
        .read(H5::PredType::NATIVE_INT, &attr_num_mutations);
    arg_archive_file.openAttribute("offset").read(H5::PredType::NATIVE_INT, &attr_offset);
    if (attr_arg_file_version == 1) {
      arg_archive_file.openAttribute("sequence_length")
        .read(H5::PredType::NATIVE_DOUBLE, &attr_sequence_length);
    }
    else if (attr_arg_file_version == 2) {
      // temporary treatment for v2 files
      arg_real_t attr_start, attr_end;
      arg_archive_file.openAttribute("start")
        .read(H5::PredType::NATIVE_DOUBLE, &attr_start);
      arg_archive_file.openAttribute("end")
        .read(H5::PredType::NATIVE_DOUBLE, &attr_end);
      attr_sequence_length = attr_end - attr_start;

    }
    else {
      std::cerr << "Unknown version of arg file" << std::endl;
    }
    
    // std::cout << "### Num nodes: " << attr_num_nodes << '\n';
    // std::cout << "### Num mutations: " << attr_num_mutations << '\n';

    std::vector<uint8_t> is_sample(attr_num_nodes);
    std::deque<bool> is_sample_deque;
    std::vector<arg_real_t> node_heights(attr_num_nodes);
    std::vector<int> edge_ids(attr_num_edges * 2);
    std::vector<arg_real_t> edge_ranges(attr_num_edges * 2);

    H5::DataSet is_sample_ds = arg_archive_file.openDataSet("flags");
    is_sample_ds.read(is_sample.data(), H5::PredType::NATIVE_HBOOL);
    is_sample_deque.assign(is_sample.begin(), is_sample.end());

    H5::DataSet times_ds = arg_archive_file.openDataSet("times");
    times_ds.read(node_heights.data(), H5::PredType::NATIVE_DOUBLE);
    // std::vector<arg_real_t> first_ten_heights(node_heights.begin()+ 10000, node_heights.begin() +
    // 10010); for(const auto& height : first_ten_heights) {
    //   std::cout << height << '\n';
    // }

    H5::DataSet edge_ids_ds = arg_archive_file.openDataSet("edge_ids");
    edge_ids_ds.read(edge_ids.data(), H5::PredType::NATIVE_INT);

    H5::DataSet edge_ranges_ds = arg_archive_file.openDataSet("edge_ranges");
    edge_ranges_ds.read(edge_ranges.data(), H5::PredType::NATIVE_DOUBLE);

    // Post process into vector of pairs as per tskit-style ARG constructor
    assert(edge_ids.size() % 2 == 0);
    std::vector<std::pair<int, int>> edge_ids_paired;
    edge_ids_paired.reserve(edge_ids.size() / 2);
    for (auto i = 0u; i < edge_ids.size(); i += 2)
      edge_ids_paired.emplace_back(edge_ids[i], edge_ids[i + 1]);

    assert(edge_ranges.size() % 2 == 0);
    std::vector<std::pair<arg_real_t, arg_real_t>> edge_ranges_paired;
    edge_ranges_paired.reserve(edge_ranges.size() / 2);
    for (auto i = 0u; i < edge_ranges.size(); i += 2)
      edge_ranges_paired.emplace_back(edge_ranges[i], edge_ranges[i + 1]);
    
    if (trim_start > 0. || trim_end < attr_sequence_length) {
      if (trim_start > trim_end) {
        throw std::logic_error("trim start is after trim end");
      }
      if (truncation_height != std::numeric_limits<arg_real_t>::max()) {
        throw std::logic_error("cannot trim and truncate at the same time, for now.");
      }
      vector<int> node_is_in_range(attr_num_nodes, 0);
      arg_real_t expected_num_edges = (arg_real_t (attr_num_edges)) * (trim_end - trim_start)/attr_sequence_length;
      vector<std::array<int, 2UL>> edge_ids_trimmed;
      edge_ids_trimmed.reserve(int (expected_num_edges));
      vector<std::array<arg_real_t, 2UL>> edge_ranges_trimmed;
      edge_ranges_trimmed.reserve(int (expected_num_edges));

      for (int edge_id = 0; edge_id < edge_ranges_paired.size(); edge_id++) {
        auto const& edge_start = edge_ranges_paired.at(edge_id).first;
        auto const& edge_end = edge_ranges_paired.at(edge_id).second;
        if (edge_start < trim_end && edge_end > trim_start) {
          std::array<arg_real_t, 2UL> edge_range{std::max(edge_start - trim_start, arg_real_t(0)), std::min(edge_end, trim_end) - trim_start};
          std::array<int, 2UL> edge_ids{edge_ids_paired.at(edge_id).first, edge_ids_paired.at(edge_id).second};
          node_is_in_range[edge_ids_paired.at(edge_id).first] = 1;
          node_is_in_range[edge_ids_paired.at(edge_id).second] = 1;
          edge_ids_trimmed.emplace_back(edge_ids);
          edge_ranges_trimmed.emplace_back(edge_range);
        }
      }
      edge_ids_trimmed.shrink_to_fit();
      edge_ranges_trimmed.shrink_to_fit();
      int num_nodes_in_range = std::count(node_is_in_range.begin(), node_is_in_range.end(), 1);
      vector<int> reassigned_node_id(attr_num_nodes, 0);
      // new node id for a node is the number of in range nodes before it -- should preserve order and
      // begin with leaf nodes
      std::partial_sum(node_is_in_range.begin(), node_is_in_range.end(), reassigned_node_id.begin());
      for (auto& elem : reassigned_node_id) {
        elem -= 1; // to make ids 0-based index
      }
      for (auto& elem : edge_ids_trimmed) {
        int new_child_id = reassigned_node_id[elem[0]];
        int new_parent_id = reassigned_node_id[elem[1]];
        std::array<int, 2UL> new_ids{new_child_id, new_parent_id};
        elem = new_ids;
      }

      // assemble all ARG data for the deserialisation constructor
      vector<arg_real_t> node_heights_trimmed(num_nodes_in_range, 0);
      vector<std::array<arg_real_t, 2UL>> node_bounds_trimmed(num_nodes_in_range);
      vector<uint8_t> is_sample_trimmed(num_nodes_in_range, false);

      for (int old_node_id = 0; old_node_id < node_is_in_range.size(); old_node_id++) {
        if (node_is_in_range[old_node_id]) {
          int new_node_id = reassigned_node_id[old_node_id];
          is_sample_trimmed[new_node_id] = uint8_t (is_sample.at(old_node_id));
          node_heights_trimmed[new_node_id] = node_heights[old_node_id];
          std::array<arg_real_t, 2UL> current_node_bound{std::max(0 - trim_start, arg_real_t(0)), std::min(attr_sequence_length, trim_end) -  trim_start};
          node_bounds_trimmed[new_node_id] = current_node_bound;
        }
      }

      node_heights.clear();
      is_sample.clear();
      edge_ids.clear();
      edge_ranges.clear();

      DeserializationParams dp{};
      dp.chromosome = attr_chromosome;
      dp.start = std::max(0 - trim_start, arg_real_t(0));
      dp.end = std::min(trim_end, attr_sequence_length) - trim_start;
      dp.offset = attr_offset + trim_start;
      dp.threaded_samples = num_nodes_in_range;
      dp.reserved_samples = -1;
      dp.num_nodes = num_nodes_in_range;
      ARG trimmed_arg = ARG(dp);
      trimmed_arg.deserialize_add_nodes(node_heights_trimmed, is_sample_trimmed, node_bounds_trimmed);
      trimmed_arg.deserialize_add_edges(edge_ids_trimmed, edge_ranges_trimmed);
      return trimmed_arg;

    }

    // std::vector<std::pair<int, int>> first_ten_ids(edge_ids_paired.begin()+ 10000,
    // edge_ids_paired.begin() + 10010); for(const auto& id : first_ten_ids) {
    //   std::cout << id.first << "," << id.second << '\n';
    // }
    // std::vector<std::pair<arg_real_t, arg_real_t>> first_ten_ranges(edge_ranges_paired.begin()+
    // 10000, edge_ranges_paired.begin() + 10010); for(const auto& range : first_ten_ranges) {
    //   std::cout << range.first << ", " << range.second <<  '\n';
    // }
    // begin truncation rewrite
    ARG_data arg_data;
    arg_data.node_heights = node_heights;
    arg_data.is_sample_deque = is_sample_deque;
    arg_data.edge_ids_paired = edge_ids_paired;
    arg_data.edge_ranges_paired = edge_ranges_paired;

    arg_data = truncate_arg_data(arg_data, truncation_height);

    ARG arg(0, attr_sequence_length, arg_data.node_heights, arg_data.is_sample_deque,
            arg_data.edge_ids_paired, arg_data.edge_ranges_paired);
    arg.set_offset(attr_offset);
    arg.set_chromosome(attr_chromosome);
    return arg;
  } // end of try block

  // catch failure caused by the H5File operations
  catch (const H5::FileIException& error) {
    H5::FileIException::printErrorStack();
  }

  // catch failure caused by the DataSet operations
  catch (const H5::DataSetIException& error) {
    H5::DataSetIException::printErrorStack();
  }
  throw std::runtime_error("ARG deserialization failed.");
}

arg_utils::ARG_data arg_utils::truncate_arg_data(ARG_data& arg_data, arg_real_t truncation_height, bool split_root) {
  if (truncation_height >= std::numeric_limits<arg_real_t>::max()) {
    return arg_data;
  }

  int new_node_index = 0;
  std::deque<bool> truncated_is_sample_deque;
  std::vector<arg_real_t> truncated_node_heights;
  truncated_node_heights.reserve(arg_data.node_heights.size());
  std::unordered_map<int, int> old_to_new_node_id;
  old_to_new_node_id.reserve(arg_data.node_heights.size());
  std::vector<std::pair<int, int>> truncated_edge_ids_paired;
  std::vector<std::pair<arg_real_t, arg_real_t>> truncated_edge_ranges_paired;

  for (auto i = 0u; i < arg_data.node_heights.size(); i++) {
    if (arg_data.node_heights[i] <= truncation_height) {
      old_to_new_node_id.emplace(i, new_node_index);
      truncated_node_heights.emplace_back(arg_data.node_heights[i]);
      truncated_is_sample_deque.emplace_back(arg_data.is_sample_deque[i]);
      new_node_index++;
    }
  }
  // OLD BEHAVIOUR
  if (!split_root) {
    // one parent root node for all edges crossing the truncation cutoff
    // put the top node epsilon higher than the cutoff point, just in case we have nodes there... (perhaps unnecessary)
    truncated_node_heights.emplace_back(truncation_height + std::numeric_limits<arg_real_t>::epsilon());
    truncated_is_sample_deque.emplace_back(0);
    int top_node_index = new_node_index;

    for (auto i = 0u; i < arg_data.edge_ids_paired.size(); i++) {
      int child_node_idx = arg_data.edge_ids_paired[i].first;
      int parent_node_idx = arg_data.edge_ids_paired[i].second;
      if (arg_data.node_heights[child_node_idx] <= truncation_height) {
        truncated_edge_ranges_paired.emplace_back(arg_data.edge_ranges_paired[i]);
        truncated_edge_ids_paired.emplace_back(
            old_to_new_node_id[child_node_idx],
            arg_data.node_heights[parent_node_idx] <= truncation_height
                ? old_to_new_node_id[parent_node_idx]
                : top_node_index);
      }
    }
  }
  else {
    // multiple parent root nodes for edges crossing the truncation cutoff
    // two edges share a root node if their ranges overlap
    std::set<arg_real_t> split_edge_endpoints;
    for (auto i = 0u; i < arg_data.edge_ids_paired.size(); i++) {
      int child_node_idx = arg_data.edge_ids_paired[i].first;
      int parent_node_idx = arg_data.edge_ids_paired[i].second;
      if (arg_data.node_heights[child_node_idx] <= truncation_height){ 
        if (arg_data.node_heights[parent_node_idx] > truncation_height) {
          split_edge_endpoints.emplace(arg_data.edge_ranges_paired[i].first);
          split_edge_endpoints.emplace(arg_data.edge_ranges_paired[i].second);
        }
          else if (arg_data.node_heights[parent_node_idx] <= truncation_height) {
          truncated_edge_ranges_paired.emplace_back(arg_data.edge_ranges_paired[i]);
          truncated_edge_ids_paired.emplace_back(
            old_to_new_node_id[child_node_idx], old_to_new_node_id[parent_node_idx]);
        }
      }
    }

    // construct new root nodes
    for (int i = 0; i < split_edge_endpoints.size() - 1; i++) {
      truncated_node_heights.emplace_back(truncation_height + std::numeric_limits<arg_real_t>::epsilon());
      truncated_is_sample_deque.emplace_back(0);
    }

    // std::cout << "done inserting " << split_edge_endpoints.size() - 1 << " new root nodes" << std::endl;
    // std::cout << "root split at ";
    // for (auto split : split_edge_endpoints) {
    //   std::cout << split << ", ";
    // }
    // std::cout << std::endl;

    // connect children nodes of truncated edges to the corresponding root nodes
    // // first convert endpoints into vector so we have constant time indexing
    // std::vector<arg_real_t> split_points(split_edge_endpoints.begin(), split_edge_endpoints.end());
    for (auto i = 0u; i < arg_data.edge_ids_paired.size(); i++) {
      int child_node_idx = arg_data.edge_ids_paired[i].first;
      int parent_node_idx = arg_data.edge_ids_paired[i].second;
      if (arg_data.node_heights[child_node_idx] <= truncation_height && 
          arg_data.node_heights[parent_node_idx] > truncation_height) {
        arg_real_t edge_start = arg_data.edge_ranges_paired[i].first;
        arg_real_t edge_end = arg_data.edge_ranges_paired[i].second;
        auto edge_start_it = split_edge_endpoints.find(edge_start);
        auto edge_end_it = split_edge_endpoints.upper_bound(edge_end);
        if (edge_start_it == split_edge_endpoints.end()) {
          throw std::logic_error("edge start is not recorded as a split position??");
        }
        edge_end_it--;


        int top_node_idx_offset = std::distance(split_edge_endpoints.begin(), edge_start_it);
        int offset = 0;
        for (auto offset_it = edge_start_it; offset_it != edge_end_it; offset_it++){
          arg_real_t current_edge_start = *(offset_it);
          offset_it++;
          if (offset_it == split_edge_endpoints.end()) {break;}
          arg_real_t current_edge_end = *(offset_it);
          offset_it--;
          int current_top_node_idx = new_node_index + top_node_idx_offset + offset;
          truncated_edge_ids_paired.emplace_back(old_to_new_node_id[child_node_idx], current_top_node_idx);
          // std::cout << "adding edge from " << old_to_new_node_id[child_node_idx] << " to " << current_top_node_idx << " with range " << current_edge_start << " to " << current_edge_end << " offset " << offset << " org start, end: " << edge_start << ", " << edge_end << std::endl;
          truncated_edge_ranges_paired.emplace_back(current_edge_start, current_edge_end);
          offset++;
        }
      }
    }
  }
  ARG_data truncated;
  truncated.node_heights = truncated_node_heights;
  truncated.is_sample_deque = truncated_is_sample_deque;
  truncated.edge_ids_paired = truncated_edge_ids_paired;
  truncated.edge_ranges_paired = truncated_edge_ranges_paired;
  return truncated;
}
