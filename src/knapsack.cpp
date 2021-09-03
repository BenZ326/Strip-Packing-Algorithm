/*
 * Copyright Xiangyi Zhang 2021
 * The code may be used for academic, non-commercial purposes only.
 * Please contact me at xiangyi.zhang@polymtl.ca for questions
 * If you have improvements, please contact me!
 */
#include "knapsack.h"

#include <algorithm>
#include <iostream>

double dynamicPrg4KnapSack(const std::vector<double>& t_values,
                           const std::vector<int>& t_weights, int capacity,
                           std::vector<int>& t_selected) {
  size_t itemSize = t_values.size() + 1;
  size_t weightSize = capacity + 1;
  double** valueMatrix = new double*[weightSize];
  for (size_t i = 0; i < weightSize; ++i) {
    valueMatrix[i] = new double[itemSize];
    for (size_t j = 0; j < itemSize; ++j) {
      valueMatrix[i][j] = 0;
    }
  }
  for (size_t j = 1; j < itemSize; ++j) {
    for (size_t i = 1; i < weightSize; ++i) {
      if (t_weights[j - 1] > i) {
        valueMatrix[i][j] = valueMatrix[i][j - 1];
      } else {
        valueMatrix[i][j] = std::max(
            (valueMatrix[(i - t_weights[j - 1])][j - 1] + t_values[j - 1]),
            double(valueMatrix[i][j - 1]));
      }
    }
  }
  double result = valueMatrix[weightSize - 1][itemSize - 1];
  double tmpRes = result;
  // provide the solution
  int W = capacity;
  for (size_t j = itemSize - 1; j-- > 0 && tmpRes > 0;) {
    if (std::abs(valueMatrix[W][j] - tmpRes) < tolerance)
      continue;
    else {
      t_selected.push_back(j);
      tmpRes -= t_values[j];
      W -= t_weights[j];
    }
  }
  for (size_t i = 0; i < weightSize; ++i) {
    delete[] valueMatrix[i];
  }
  delete[] valueMatrix;
  return result;
}
