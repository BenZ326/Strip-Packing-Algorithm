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
