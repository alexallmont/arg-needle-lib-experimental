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

#ifndef ARG_NEEDLE_LIB_DESERIALIZATION_PARAMS_HPP
#define ARG_NEEDLE_LIB_DESERIALIZATION_PARAMS_HPP

#include "types.hpp"

/**
 * @brief Struct containing parameters for ARG deserialization
 *
 */
struct DeserializationParams {
  arg_real_t start{};
  arg_real_t end{};
  int num_nodes{};
  int offset{};
  int chromosome{};
  int threaded_samples{};
  int reserved_samples{};
};

#endif // ARG_NEEDLE_LIB_DESERIALIZATION_PARAMS_HPP
