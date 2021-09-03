#include "BLEU.h"

/*
processing for the y-check algorithm
*/

const std::vector<const StripPacking::item*>
StripPacking::BLEU::preprocess4yCheck(
    int& t_Width, const std::vector<const item*>& t_InterestItems,
    std::vector<coordinate>& t_Cords, const int t_Height) const {
  auto Items =
      this->preprocessedFirst4yCheck(t_InterestItems, t_Cords, t_Width);
  this->preprocessedSecond4yCheck(Items, t_Cords, t_Width);
  this->preprocessedThird4yCheck(Items, t_Cords, t_Width);
  std::vector<const item*> result;
  for (const auto& it : Items) {
    result.push_back(std::move(it));
  }
  return result;
}

const std::vector<StripPacking::item*>
StripPacking::BLEU::preprocessedFirst4yCheck(
    const std::vector<const item*>& t_InterestItems,
    std::vector<coordinate>& t_Cords, const int t_Width) const {
  std::vector<item*> allItems;
  // copy all the items
  for (const auto& it : t_InterestItems) {
    item* newItem = new item(*it);
    allItems.push_back(newItem);
    newItem->subItems.clear();
  }
  std::set<int> coveredItems;
  auto maps = this->getLeftsAndRights(allItems, t_Cords, coveredItems);
  while (true) {
    auto& leftItems = maps[0];
    auto& rightItems = maps[1];
    bool exit = true;
    for (const auto& it : leftItems) exit &= (it.second.empty());
    for (const auto& it : rightItems) exit &= (it.second.empty());
    if (exit) break;
    std::sort(allItems.begin(), allItems.end(), compareItemByxCords(t_Cords));
    for (size_t i = 0; i < allItems.size(); ++i) {
      item* tmpItem = allItems[i];
      if (coveredItems.find(tmpItem->idx) != coveredItems.end()) continue;
      auto& lefts = leftItems.find(tmpItem->idx)->second;
      auto& rights = rightItems.find(tmpItem->idx)->second;
      // merge left
      bool leftMergeable =
          mergeItems4yCheck(allItems, tmpItem, leftItems, rightItems, lefts,
                            t_Cords, true, t_Width);
      // merge right
      bool rightMergeable =
          mergeItems4yCheck(allItems, tmpItem, leftItems, rightItems, rights,
                            t_Cords, false, t_Width);
      if (leftMergeable || rightMergeable) {
        maps = this->getLeftsAndRights(allItems, t_Cords, coveredItems);
        break;
      }
    }
  }
  return allItems;
}

/*
Lift Item width
Args:
        t_allItems: all the items after the first processed technique
        t_Cords: the coordinates of each items in t_allItems
*/
void StripPacking::BLEU::preprocessedSecond4yCheck(
    std::vector<item*>& t_allItems, std::vector<coordinate>& t_Cords,
    const int t_binWidth) const {
  auto maps = this->getLeftsAndRights(t_allItems, t_Cords, std::set<int>());
  std::sort(t_allItems.begin(), t_allItems.end(), compareItemByWidthLess);
  for (const auto& it : t_allItems) {
    auto leftItems = maps[0].find(it->idx)->second;
    auto rightItems = maps[1].find(it->idx)->second;
    auto l_j = 0;
    for (const auto& jt : leftItems) {
      if (l_j < t_Cords[jt->idxHelper].x + jt->width) {
        l_j = t_Cords[jt->idxHelper].x + jt->width;
      }
    }
    auto r_j = t_binWidth;
    for (const auto& jt : rightItems) {
      if (r_j > t_Cords[jt->idxHelper].x) {
        r_j = t_Cords[jt->idxHelper].x;
      }
    }
    t_Cords[it->idxHelper].x = l_j;
    it->width = r_j - l_j;
  }
}

/*
Lift Item width
Args:
        t_allItems: all the items after the second processed technique
        t_Cords: the coordinates of each items in t_allItems
        t_binWidth : the width of the bin
*/
void StripPacking::BLEU::preprocessedThird4yCheck(
    std::vector<item*>& t_allItems, std::vector<coordinate>& t_Cords,
    int& t_binWidth) const {
  std::map<int, item*> idxHelper_Items;
  for (const auto& it : t_allItems) {
    if (t_Cords[it->idxHelper].x != -1)  // dead item
    {
      idxHelper_Items.insert(std::pair<int, item*>(it->idxHelper, it));
    }
  }
  int cur = 0;
  int probe = 1;
  auto mergeCol = [](const int t_sCol, const int t_eCol,
                     std::vector<item*>& t_ItemsInCol,
                     std::vector<item*>& t_ItemsAfterEndCol,
                     std::vector<coordinate>& t_Cords, int& t_binWidth) {
    // all items in ItemsInColshould be reduced by t_eCol - t_sCol in width
    // x coordinate of the items in ItemsAfterEndCol should should be reduced by
    // t_eCol - t_sCol t_binWidth -= t_eCol - t_sCol
    auto diff = t_eCol - t_sCol;
    for (auto& it : t_ItemsInCol) it->width -= diff;
    for (const auto& it : t_ItemsAfterEndCol) t_Cords[it->idxHelper].x -= diff;
    t_binWidth -= diff;
    return;
  };
  while (cur < t_binWidth) {
    const auto&& occItemsCur = this->getItemsByCol(cur, t_allItems, t_Cords);
    // locate the right-most columns with the same items as current column
    while (probe < t_binWidth) {
      // yes, go on probe
      // no, merge columns between [cur, probe-1]
      const auto&& occItemsProb =
          this->getItemsByCol(probe, t_allItems, t_Cords);
      // if current column can be merged with the column being probed, continue,
      // go on probing
      if (occItemsCur == occItemsProb) {
        probe++;
      } else
        break;
    }
    // found the right-most column
    if (probe - cur > 1) {
      // merge [cur, probe - 1]
      std::vector<item*> ItemsInCol;
      for (const auto& it : occItemsCur)
        ItemsInCol.push_back(idxHelper_Items.find(it)->second);
      std::vector<item*> ItemsAfterEndCol;
      for (const auto& it : t_allItems) {
        if (t_Cords[it->idxHelper].x >= probe) ItemsAfterEndCol.push_back(it);
      }
      mergeCol(cur, probe - 1, ItemsInCol, ItemsAfterEndCol, t_Cords,
               t_binWidth);
      cur++;
      probe = cur + 1;
      continue;
    }
    // update cur and probe
    cur = probe;
    probe++;
  }
  // assign the new coordinates
  std::vector<coordinate> newCords(t_allItems.size(), coordinate(0, 0));
  int idxHelper = 0;
  for (const auto& it : t_allItems) {
    newCords[idxHelper] = t_Cords[it->idxHelper];
    it->idxHelper = idxHelper;
    idxHelper++;
  }
  t_Cords = std::move(newCords);
}

/*
The following is well described in the 3.2 section of the paper
1) first step:
        if leftItems is empty, exit function with false returned
        check all the left items to see if h_i <= h_j for any i
                if yes, try to merge the items, if mergeable, then return true,
if not go to the second step 2) second step: get t_startColumn, get the items
that are packed on t_startColumn with the largest w. then only remains item i in
the left that satisfy: p_i >= t_startColumn + w, p_i + w_i <= p_j then go to the
first step
*/
bool StripPacking::BLEU::mergeItems4yCheck(
    std::vector<item*>& t_allItems, item* t_i,
    std::map<int, std::list<item*>>& t_leftItems,
    std::map<int, std::list<item*>>& t_rightItems, std::list<item*>& t_Items,
    std::vector<coordinate>& t_Cords, const bool t_left,
    const int t_binWidth) const {
  if (t_Items.empty()) return false;
  // check lefts
  auto& p_js = t_Cords[t_i->idxHelper].x;
  if (t_left) {
    int startColumn = 0;
    int maxWidth = 1;
    while (!t_Items.empty()) {
      // step 1
      std::vector<int> allHeights;
      std::vector<const item*> transferredItems;
      std::vector<coordinate> transferredCords;
      for (const auto& it : t_Items) allHeights.push_back(it->height);
      if (*(std::max_element(allHeights.begin(), allHeights.end())) <=
          t_i->height) {
        this->transferItemsAndCords4YEnumeration(
            t_Items, t_Cords, transferredItems, transferredCords);
        if (this->yCheckEnumerationTree(
                transferredItems, transferredCords, t_i->height,
                t_Cords[t_i->idxHelper].x - (startColumn + maxWidth - 1)) ==
            solutionStatus::feasible) {
          this->merging(t_i, t_allItems, t_Cords, startColumn, maxWidth,
                        t_Items, true);
          this->releaseTmpItems(transferredItems);
          return true;
        } else
          this->releaseTmpItems(transferredItems);
      }
      startColumn = this->getFirstColumn(t_Items, t_Cords);
      maxWidth = this->getMaxWidth(startColumn, t_Items, t_Cords, true);
      while (!this->checkSeparable(t_Items, t_Cords, startColumn, maxWidth,
                                   true))  // means it can be spearated
      {
        startColumn = this->getFirstColumn(startColumn, t_Items, t_Cords);
        maxWidth = this->getMaxWidth(startColumn, t_Items, t_Cords, true);
        if (startColumn < 0) return false;  // when t_Items is empty
      }
      // separate the items in the left side
      auto ptr = t_Items.begin();
      while (ptr != t_Items.end()) {
        auto& p_is = t_Cords[(*ptr)->idxHelper].x;
        if (!(startColumn + maxWidth < p_is + (*ptr)->width &&
              p_is + (*ptr)->width <= p_js)) {
          ptr = t_Items.erase(ptr);
        } else
          ptr++;
      }
    }
    return false;
  }
  // check rights
  else {
    int LastColumn = t_binWidth - 1;
    int maxWidth = 0;
    while (!t_Items.empty()) {
      // step 1
      std::vector<int> allHeights;
      std::vector<const item*> transferredItems;
      std::vector<coordinate> transferredCords;
      for (const auto& it : t_Items) allHeights.push_back(it->height);
      if (*(std::max_element(allHeights.begin(), allHeights.end())) <=
          t_i->height) {
        this->transferItemsAndCords4YEnumeration(
            t_Items, t_Cords, transferredItems, transferredCords);
        if (this->yCheckEnumerationTree(
                transferredItems, transferredCords, t_i->height,
                LastColumn - maxWidth -
                    (t_Cords[t_i->idxHelper].x + t_i->width - 1)) ==
            solutionStatus::feasible) {
          this->merging(t_i, t_allItems, t_Cords, LastColumn, maxWidth, t_Items,
                        false);
          this->releaseTmpItems(transferredItems);
          return true;
        } else
          this->releaseTmpItems(transferredItems);
      }
      LastColumn = this->getLastColumn(t_Items, t_Cords);
      maxWidth = this->getMaxWidth(LastColumn, t_Items, t_Cords, false);
      while (!this->checkSeparable(t_Items, t_Cords, LastColumn, maxWidth,
                                   false))  // means it can be spearated
      {
        LastColumn = this->getLastColumn(LastColumn, t_Items, t_Cords);
        maxWidth = this->getMaxWidth(LastColumn, t_Items, t_Cords, false);
        if (LastColumn < 0) return false;  // when t_Items is empty
      }
      // separate the items in the left side
      auto ptr = t_Items.begin();
      while (ptr != t_Items.end()) {
        auto& p_is = t_Cords[(*ptr)->idxHelper].x;
        if (!(p_is + (*ptr)->width <= LastColumn - maxWidth && p_is < p_js)) {
          ptr = t_Items.erase(ptr);
        } else
          ptr++;
      }
    }
    return false;
  }
}

const bool StripPacking::BLEU::checkSeparable(
    const std::list<item*>& t_Items, const std::vector<coordinate>& t_Cords,
    const int t_startColumn, const int t_maxWidth, const bool t_left) const {
  if (t_left) {
    for (const auto it : t_Items) {
      auto& p_is = t_Cords[it->idxHelper].x;
      if (p_is < t_startColumn + t_maxWidth &&
          p_is + it->width > t_startColumn + t_maxWidth)
        return false;
    }
    return true;
  } else {
    for (const auto it : t_Items) {
      auto& p_is = t_Cords[it->idxHelper].x;
      if (p_is + it->width > t_startColumn - t_maxWidth &&
          p_is < t_startColumn - t_maxWidth)
        return false;
    }
    return true;
  }
}

const int StripPacking::BLEU::getFirstColumn(
    const std::list<item*>& t_lefts,
    const std::vector<coordinate>& t_Cords) const {
  if (t_lefts.empty()) return -1;
  std::vector<int> xCords;
  for (const auto& it : t_lefts) {
    xCords.push_back(t_Cords[it->idxHelper].x);
  }
  return xCords.empty() ? -1 : *std::min_element(xCords.begin(), xCords.end());
}

const int StripPacking::BLEU::getFirstColumn(
    const int t_prevStartColumn, const std::list<item*>& t_lefts,
    const std::vector<coordinate>& t_Cords) const {
  if (t_lefts.empty()) return -1;
  std::vector<int> xCords;
  for (const auto& it : t_lefts) {
    if (t_Cords[it->idxHelper].x > t_prevStartColumn)
      xCords.push_back(t_Cords[it->idxHelper].x);
  }
  return (xCords.empty() ? -1
                         : *std::min_element(xCords.begin(), xCords.end()));
}

const int StripPacking::BLEU::getLastColumn(
    const std::list<item*>& t_rights,
    const std::vector<coordinate>& t_Cords) const {
  if (t_rights.empty()) return -1;
  std::vector<int> xCords;
  for (const auto& it : t_rights) {
    xCords.push_back(t_Cords[it->idxHelper].x + it->width - 1);
  }
  return xCords.empty() ? -1 : *std::max_element(xCords.begin(), xCords.end());
}

const int StripPacking::BLEU::getLastColumn(
    const int t_prevLastColumn, const std::list<item*>& t_rights,
    const std::vector<coordinate>& t_Cords) const {
  if (t_rights.empty()) return -1;
  std::vector<int> xCords;
  for (const auto& it : t_rights) {
    if (t_Cords[it->idxHelper].x + it->width - 1 < t_prevLastColumn)
      xCords.push_back(t_Cords[it->idxHelper].x + it->width - 1);
  }
  return xCords.empty() ? -1 : *std::max_element(xCords.begin(), xCords.end());
}

/*
Given a column, t_column, find the item that has the largest width of which the
first peice being placed at column t_column
*/
const int StripPacking::BLEU::getMaxWidth(
    const int t_column, const std::list<item*>& t_Items,
    const std::vector<coordinate>& t_Cords, bool t_left) const {
  if (t_left) {
    std::vector<int> widths;
    for (const auto& it : t_Items) {
      if (t_Cords[it->idxHelper].x == t_column) widths.push_back(it->width);
    }
    return *std::max_element(widths.begin(), widths.end());
  } else {
    std::vector<int> widths;
    for (const auto& it : t_Items) {
      if (t_Cords[it->idxHelper].x + it->width - 1 == t_column)
        widths.push_back(it->width);
    }
    return *std::max_element(widths.begin(), widths.end());
  }
}
void StripPacking::BLEU::transferItemsAndCords4YEnumeration(
    const std::list<item*>& t_OrigItems,
    const std::vector<coordinate>& t_OrigCords,
    std::vector<const item*>& t_Items, std::vector<coordinate>& t_Cords) const {
  int idxHelper = 0;
  std::vector<int> oldCords;
  for (const auto& it : t_OrigItems) {
    coordinate cord(t_OrigCords[it->idxHelper].x, t_OrigCords[it->idxHelper].y);
    t_Cords.push_back(cord);
    oldCords.push_back(cord.x);
    const item* tmp = new item(it->idx, it->width, it->height, idxHelper++);
    t_Items.push_back(tmp);
  }
  // transfer old coorindate to new coordinate
  std::sort(oldCords.begin(), oldCords.end());
  std::map<int, int> oldNewCords;
  int idx = 0;
  for (const auto& it : oldCords) {
    if (oldNewCords.find(it) == oldNewCords.end()) {
      oldNewCords.insert(std::pair<int, int>(it, idx++));
    }
  }
  for (size_t i = 0; i < t_Cords.size(); ++i) {
    t_Cords[i].x = oldNewCords.find(t_Cords[i].x)->second;
  }
}

void StripPacking::BLEU::merging(item* t_i, std::vector<item*>& t_allItems,
                                 std::vector<coordinate>& t_Cords,
                                 const int t_startColumn, const int t_maxWidth,
                                 std::list<item*>& t_Items,
                                 const bool t_left) const {
  if (t_left) {
    // merge
    // 1) update t_i
    // input, t_items, t_Cords, t_allItems, startColumn, maxWidth,
    for (const auto& it : t_Items) {
      t_i->subItems.push_back(it);
      t_Cords[it->idxHelper].x = -1;
    }
    t_i->width += (t_Cords[t_i->idxHelper].x - t_startColumn - t_maxWidth);
    t_Cords[t_i->idxHelper].x = t_startColumn + t_maxWidth;
    auto ptr = t_allItems.begin();
    // 2) update allItems
    while (ptr != t_allItems.end()) {
      bool found = false;
      for (const auto& it : t_i->subItems) {
        if ((*ptr)->idx == it->idx) {
          found = true;
          ptr = t_allItems.erase(ptr);
          break;
        }
      }
      if (!found) ptr++;
    }
    // 3) update maps
    t_Items.clear();
  } else {
    // merge
    // 1) update t_i
    // input, t_items, t_Cords, t_allItems, startColumn, maxWidth,
    for (const auto& it : t_Items) {
      t_i->subItems.push_back(it);
      t_Cords[it->idxHelper].x = -1;
    }
    t_i->width += (t_startColumn - t_maxWidth -
                   (t_Cords[t_i->idxHelper].x + t_i->width - 1));
    auto ptr = t_allItems.begin();
    // 2) update allItems
    while (ptr != t_allItems.end()) {
      bool found = false;
      for (const auto& it : t_i->subItems) {
        if ((*ptr)->idx == it->idx) {
          found = true;
          ptr = t_allItems.erase(ptr);
          break;
        }
      }
      if (!found) ptr++;
    }
    // 3) update maps
    t_Items.clear();
  }
}

void StripPacking::BLEU::mergeItems(item* t_i, std::list<item*>& t_Items,
                                    std::map<int, std::list<item*>>& t_allItems,
                                    std::vector<coordinate>& t_Cords) const {}

/*
if an item is in t_coveredItems, then it has no left items and right items
*/
const std::vector<std::map<int, std::list<StripPacking::item*>>>
StripPacking::BLEU::getLeftsAndRights(
    const std::vector<item*>& t_allItems,
    const std::vector<coordinate>& t_Cords,
    const std::set<int>& t_coveredItems) const {
  std::vector<std::map<int, std::list<item*>>> res;
  std::map<int, std::list<item*>> leftItems;
  std::map<int, std::list<item*>> rightItems;
  for (const auto& j_it : t_allItems) {
    leftItems.insert(
        std::pair<int, std::list<item*>>(j_it->idx, std::list<item*>()));
    rightItems.insert(
        std::pair<int, std::list<item*>>(j_it->idx, std::list<item*>()));
    if (t_coveredItems.find(j_it->idx) != t_coveredItems.end()) continue;
    for (const auto& i_it : t_allItems) {
      // left items
      if (j_it->idx == i_it->idx) continue;

      auto p_js = t_Cords[j_it->idxHelper].x;
      auto p_is = t_Cords[i_it->idxHelper].x;
      if (p_is + i_it->width <= p_js) {
        leftItems.find(j_it->idx)->second.push_back(i_it);
      }
      // right items
      if (p_js + j_it->width <= p_is) {
        rightItems.find(j_it->idx)->second.push_back(i_it);
      }
    }
  }
  res.push_back(leftItems);
  res.push_back(rightItems);
  return res;
}

void StripPacking::BLEU::releaseTmpItems(
    std::vector<const item*>& t_Items) const {
  for (const auto& it : t_Items) {
    if (!it->subItems.empty()) {
      this->releaseTmpItems(const_cast<item*>(it)->subItems);
    }
  }
}
