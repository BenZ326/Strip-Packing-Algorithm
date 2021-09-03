/*
 * Copyright Xiangyi Zhang 2021
 * The code may be used for academic, non-commercial purposes only.
 * Please contact me at xiangyi.zhang@polymtl.ca for questions
 * If you have improvements, please contact me!
 */
#pragma once

#include "spp.h"
namespace StripPacking {
enum skylineSelectionMode { leftBottom, bestFit };
struct Skyline {
 public:
  Skyline(const int t_stripW)
      : corX(0), corY(0), length(t_stripW), next(nullptr), prev(nullptr) {}
  Skyline() {}
  Skyline(const bool t_dummy)
      : corX(-1), corY(999999), length(0), next(nullptr), prev(nullptr) {}
  int corX;
  int corY;
  int length;
  Skyline* next;  // from left to right
  Skyline* prev;
};
// select a skyline for an item to place according to the given mode
Skyline* selectSkyline(const Skyline* t_head,
                       const skylineSelectionMode& t_mode);
// add an item over the selected skyline, see the relaxationMode
void addItemOverSkyline(Skyline* t_skyline, const item* t_item);

void detectAndMergeSkylines(Skyline* t_skyline);
void removeSkyline(Skyline* t_skyline);
void liftSkyline(Skyline* t_skyline);  // lift a skyline to the adjacent skyline
                                       // that has lower height

}  // namespace StripPacking
