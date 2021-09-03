/*
 * Copyright Xiangyi Zhang 2021
 * The code may be used for academic, non-commercial purposes only.
 * Please contact me at xiangyi.zhang@polymtl.ca for questions
 * If you have improvements, please contact me!
 */
#pragma once
#include <algorithm>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>
namespace StripPacking {

enum solutionStatus {
  feasible,
  infeasible,
  pending,  // time limit or node limit makes a solution pending
  numberStatus
};

enum algorithmStatus { approximate, exact, numberAlgStatus };
/*
All fast utility function and basic structure of strip packing problem
*/

/*
pieces for parallel scheduling problem with contiguity constraints
*/
class itemPiece {
 public:
  itemPiece(const std::string t_id, const int t_height)
      : id(t_id), height(t_height) {}
  const std::string id;
  const int height;
};

/*
pieces for 1CBP
*/
class itemPieceWidth {
 public:
  itemPieceWidth(const int t_id, const int t_width)
      : id(t_id), width(t_width) {}
  const int id;
  const int width;
};

struct whPair {
  whPair() {}
  whPair(int t_h, int t_w) : h(t_h), w(t_w) {}
  whPair(int t_h, int t_w, int t_x) : h(t_h), w(t_w){};
  whPair(const whPair& t_wh) : h(t_wh.h), w(t_wh.w) {}
  bool operator==(const whPair& o1) const {
    return ((this->h == o1.h) && (this->w == o1.w));
  }
  bool operator<(const whPair& wh) const {
    return (h < wh.h) || ((h == wh.h) && (w < wh.w));
  }
  whPair operator=(const whPair& wh) {
    this->h = wh.h;
    this->w = wh.w;
    // this->xCord = wh.xCord;
    return *this;
  }
  int h;
  int w;
  // int xCord;		// x coordinate used for y check
};

class item {
 public:
  item(const int t_idx, const int t_width, const int t_height);
  item(const int t_idx, const int t_width, const int t_height,
       const int t_idxHelper);
  int idx;
  int width;
  int height;
  int idxHelper;  // an alternative identifier
  std::vector<const item*> subItems;
  static std::ostringstream ss;
  bool operator<(const item& t_item) const { return this->idx < t_item.idx; }
};

struct coordinate {
 public:
  coordinate(int t_x, int t_y) : x(t_x), y(t_y){};
  bool operator==(const coordinate& cor1) const {
    return ((x == cor1.x) && (y == cor1.y));
  }
  bool operator<(const coordinate& cor1) const {
    return (x < cor1.x || (x == cor1.x) && (y < cor1.y));
  }
  coordinate& operator=(const coordinate& cor1) {
    x = cor1.x;
    y = cor1.y;
    return *this;
  }
  int x;
  int y;
};

constexpr int BigNumber = 999999;
std::set<int> computeFX(const int t_x, const int t_idx,
                        const std::vector<const item*>& t_items, bool flag);
std::set<int> computeFX(const int t_x, const int t_idx,
                        const std::vector<StripPacking::item*>& t_items,
                        bool flag);
int getMaximalHeight(const std::vector<const item*>& t_items);
/*
Instance of the SPP
*/
double solve(const std::vector<const item*>& t_allItems,
             const std::map<int, std::set<int>>& t_mapPosWidth,
             const std::map<int, std::set<int>>& t_mapPosHeight,
             const bool t_Integer);

inline const std::string getVarName(const int t_itemIdx, const int t_xPos);

const std::string getVarName(const int t_itemIdx, const int t_xPos) {
  std::string str =
      "item" + std::to_string(t_itemIdx) + "assign" + std::to_string(t_xPos);
  return str;
}

/*
Utilities for strip packing algorithms
*/
inline bool compareItemByWidth(const item* t_i, const item* t_j);
bool compareItemByWidth(const item* t_i, const item* t_j) {
  return (t_i->width > t_j->width ||
          (t_i->width == t_j->width && t_i->height > t_j->height) ||
          (t_i->width == t_j->width && t_i->height == t_j->height &&
           t_i->idx > t_j->idx));
}

inline bool compareItemByWidthLess(const item* t_i, const item* t_j);
bool compareItemByWidthLess(const item* t_i, const item* t_j) {
  return (t_i->width < t_j->width ||
          (t_i->width == t_j->width && t_i->height < t_j->height) ||
          (t_i->width == t_j->width && t_i->height == t_j->height &&
           t_i->idx < t_j->idx));
}

inline bool compareItemByHeight(const item* t_i, const item* t_j);
bool compareItemByHeight(const item* t_i, const item* t_j) {
  return t_i->height < t_j->height ||
         (t_i->height == t_j->height && t_i->idx > t_j->idx);
}

inline bool compareItemByIdx(const item* t_i, const item* t_j);
bool compareItemByIdx(const item* t_i, const item* t_j) {
  return t_i->idx < t_j->idx;
}

inline bool compareItemByArea(const item* t_i, const item* t_j);
bool compareItemByArea(const item* t_i, const item* t_j) {
  return t_i->height * t_i->width < t_j->height * t_j->width;
}

class compareItemByWHDifference {
 public:
  const int _H;
  const int _W;
  compareItemByWHDifference(const int t_H, const int t_W) : _H(t_H), _W(t_W) {}
  bool operator()(const item* t_i, const item* t_j) const {
    return std::min(_W - t_i->width, _H - t_i->height) >
           std::min(_W - t_j->width, _H - t_j->height);
  }
};

class compareItemByHeight {
  // for heap
 public:
  bool operator()(const item* t_i, const item* t_j) const {
    return t_i->height > t_j->height ||
           (t_i->height == t_j->height && t_i->idx > t_j->idx);
  }
};

class compareItemByxCords {
 public:
  std::vector<coordinate> xCords;
  compareItemByxCords(const std::vector<coordinate>& t_xCords)
      : xCords(t_xCords) {}
  bool operator()(const item* t_i, const item* t_j) const {
    return xCords[t_i->idxHelper].x > xCords[t_j->idxHelper].x ||
           (xCords[t_i->idxHelper].x == xCords[t_j->idxHelper].x &&
            t_i->idx > t_j->idx);
  }
};

int subSetSum(const std::vector<int>& t_v, const int t_limit);
/*
Args:
        Given a vector of integers and a limit. Find the subset of the integers
of which the sum is the largest but not exceeding the limit This is achieved by
a simple dynamic programming algorithm Returns: The best sum.
*/

}  // namespace StripPacking
