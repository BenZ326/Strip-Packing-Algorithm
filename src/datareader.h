/*
 * Copyright Xiangyi Zhang 2021
 * The code may be used for academic, non-commercial purposes only.
 * Please contact me at xiangyi.zhang@polymtl.ca for questions
 * If you have improvements, please contact me!
 */
#pragma once
#include <fstream>
#include <list>
#include <sstream>
#include <vector>

#include "spp.h"

/*
cut an item into pieces each of which has 1 unit of width and the same height as
the item
*/

int readData(const std::string& t_file,
             std::vector<const StripPacking::item*>& t_items);

std::vector<int> parseLine(std::string t_line);
int convertVec(std::vector<int>& tmp);
