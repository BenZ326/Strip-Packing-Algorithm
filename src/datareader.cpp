#include "datareader.h"
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include<math.h>
#include "spp.h"



int readData(const std::string& t_file, std::vector<const StripPacking::item*>& t_items)
{
	int maxWidth;
	std::ifstream ff(t_file);
	std::string line;
	if (ff.is_open())
	{
		int idx = 1;
		std::getline(ff, line);
		int n = std::stoi(line);
		while (std::getline(ff, line))
		{
			if (idx == 1)
			{
				maxWidth = std::stoi(line);
				++idx;
				continue;
			}
			auto res = parseLine(line);
			StripPacking::item* rec = new StripPacking::item(res[0], res[1], res[2]);
			t_items.push_back(rec);
			if (res[0] == n) break;
		}

	}
	else 		std::cout << "cann't open the file" << t_file << std::endl;
	return maxWidth;
}

int convertVec(std::vector<int>& tmp)
{
	int sum = 0;
	for (size_t i = 0; i < tmp.size(); ++i) sum += tmp[i] * std::pow(10, (tmp.size() - i - 1));
	tmp.clear();
	return sum;
}

std::vector<int> parseLine(std::string t_line)
{
	std::vector<int> res;
	std::string sep = " ";
	size_t pos = 0;
	std::vector<int> tmp;
	for (auto it : t_line)
	{
		if (isspace(it))
		{
			if (!tmp.empty()) res.push_back(convertVec(tmp));
			continue;
		}
		else tmp.push_back(it - '0');
	}
	if (!tmp.empty()) res.push_back(convertVec(tmp));
	return res;
}
