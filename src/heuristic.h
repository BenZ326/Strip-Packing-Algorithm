/*
 * Copyright Xiangyi Zhang 2021
 * The code may be used for academic, non-commercial purposes only.
 * Please contact me at xiangyi.zhang@polymtl.ca for questions
 * If you have improvements, please contact me!
 */
#pragma once
#include "skyline.h"
#include "spp.h"

namespace StripPacking {
class Heuristic {
 public:  // static member functions and variables
  static std::vector<coordinate>
      solutions;  // the order is the same as the input order

 public:
  const int leftBottomHeuristic(
      const std::vector<const StripPacking::item*>& t_allItems,
      const int t_binWidth);
  const int bestFitHeuristic(std::vector<const StripPacking::item*>& t_allItems,
                             const int t_binWidth);
  const bool generalBestFitHeurisitic(
      std::vector<const StripPacking::item*>& t_allItems,
      const std::vector<const StripPacking::item*>& t_Bins);
  void dumpSolution(const std::vector<const StripPacking::item*>& t_allItems);
  const int iteratedGreedy(std::vector<const StripPacking::item*>& t_allItems,
                           const int t_binWidth);

 protected:
  const StripPacking::item* findBestItem(
      std::vector<const StripPacking::item*>& t_allItems,
      const StripPacking::Skyline* t_skyline);
  const StripPacking::item* findBestItem(
      std::vector<const StripPacking::item*>& t_allItems,
      const StripPacking::Skyline* t_skyline, const StripPacking::item* t_bin);
  const int parseSol(const Skyline* t_skyline);
};

inline void insertItem(std::vector<const StripPacking::item*>& v,
                       const int t_fromPos, const int t_toPos) {
  auto obj = v[t_fromPos];
  v.erase(v.begin() + t_fromPos);
  v.insert(v.begin() + t_toPos, obj);
}
}  // namespace StripPacking
