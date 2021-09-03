#include "spp.h"
#include <algorithm>
#include <list>
#include <set>
#include <ilcplex/ilocplex.h>

/*
The paper section 2.1 (1) and (2)
t_items: all the items
t_idx: the item excluded (idx in the t_items)
flag = false, --> height 
flag = true, -->  width
*/


std::ostringstream StripPacking::item::ss;
StripPacking::item::item(const int t_idx, const int t_width, const int t_height) :idx(t_idx), width(t_width), height(t_height)
{
	subItems.clear();
}

StripPacking::item::item(const int t_idx, const int t_width, const int t_height, const int t_idxHelper) : idx(t_idx), width(t_width), height(t_height),
idxHelper(t_idxHelper)
{
	subItems.clear();
}

std::set<int> StripPacking::computeFX(const int t_x,  const int t_idx,
	const std::vector<const StripPacking::item*>& t_items, bool flag)
{
	std::vector<std::vector<int>> res;
	for (size_t i = 0; i < t_items.size(); ++i)
	{
		std::vector<int> tmp(t_x+1, 0);
		res.push_back(tmp);
	}
	for (size_t j = 1; j <= t_x; ++j) res[0][j] = 0;
	for (size_t i = 0; i < t_items.size(); ++i) res[i][0] = 1;
	int itemIdx;
	for (size_t j = 1; j <= t_x; ++j)		// W 
	{
		for (size_t i = 1; i < t_items.size(); ++i)
		{
			if (i - 1 < t_idx) itemIdx = i - 1;
			else itemIdx = i;
			if (flag)
			{
				int diff = j - t_items[itemIdx]->width;
				if ((diff) < 0) res[i][j] = res[i - 1][j];
				else res[i][j] = std::max(res[i - 1][j], res[i - 1][j - t_items[itemIdx]->width]);
			}
			else
			{
				int diff = j - t_items[itemIdx]->height;
				if ((diff) < 0) res[i][j] = res[i - 1][j];
				else res[i][j] = std::max(res[i - 1][j], res[i - 1][j - t_items[itemIdx]->height]);
			}
		}
	}
	std::set<int> possiblePositions;
	for (size_t j = 0; j <= t_x; ++j)
	{
		if (res[res.size() - 1][j]  == 1) possiblePositions.insert(j);
	}
	return possiblePositions;
}


std::set<int> StripPacking::computeFX(const int t_x, const int t_idx,
	const std::vector<StripPacking::item*>& t_items, bool flag)
{
	std::vector<std::vector<int>> res;
	for (size_t i = 0; i < t_items.size(); ++i)
	{
		std::vector<int> tmp(t_x + 1, 0);
		res.push_back(tmp);
	}
	for (size_t j = 1; j <= t_x; ++j) res[0][j] = 0;
	for (size_t i = 0; i < t_items.size(); ++i) res[i][0] = 1;
	int itemIdx;
	for (size_t j = 1; j <= t_x; ++j)		// W 
	{
		for (size_t i = 1; i < t_items.size(); ++i)
		{
			if (i - 1 < t_idx) itemIdx = i - 1;
			else itemIdx = i;
			if (flag)
			{
				int diff = j - t_items[itemIdx]->width;
				if ((diff) < 0) res[i][j] = res[i - 1][j];
				else res[i][j] = std::max(res[i - 1][j], res[i - 1][j - t_items[itemIdx]->width]);
			}
			else
			{
				int diff = j - t_items[itemIdx]->height;
				if ((diff) < 0) res[i][j] = res[i - 1][j];
				else res[i][j] = std::max(res[i - 1][j], res[i - 1][j - t_items[itemIdx]->height]);
			}
		}
	}
	std::set<int> possiblePositions;
	for (size_t j = 0; j <= t_x; ++j)
	{
		if (res[res.size() - 1][j] == 1) possiblePositions.insert(j);
	}
	return possiblePositions;
}


int StripPacking::getMaximalHeight(const std::vector<const StripPacking::item*>& t_items)
{
	int res = -1;
	for (size_t i = 0; i < t_items.size(); ++i)
	{
		res = std::max(t_items[i]->height, res);
	}
	return res;
}


/*
Build the contiguity parallel machine scheduling problem as a lower bound for the spp
*/
double StripPacking::solve(const std::vector<const StripPacking::item*>& t_allItems, const std::map<int, std::set<int>>& t_mapPosWidth,
	const std::map<int, std::set<int>>& t_mapPosHeight, const bool t_Integer)
{
	// data preparation
	std::set<int> allPositions;
	for (const auto& it : t_mapPosWidth)
		for (const auto& it2 : it.second)
			allPositions.insert(it2);
	IloEnv env;
	IloModel model(env);
	std::map<std::string, IloNumVar> allVars;
	// first constraints set
	for (const auto& it : t_allItems)
	{
		IloExpr expr(env);
		for (const auto& it2 : t_mapPosWidth.find(it->idx)->second)
		{
			auto varName = StripPacking::getVarName(it->idx, it2);
			if (t_Integer)
			{
				IloNumVar var(env, 0, 1, ILOINT, varName.c_str());
				allVars.insert(std::pair<std::string, IloNumVar>(varName, var));
				expr += var;
			}
			else
			{
				IloNumVar var(env, 0, 1, ILOFLOAT, varName.c_str());
				allVars.insert(std::pair<std::string, IloNumVar>(varName, var));
				expr += var;
			}

		}
		model.add(expr == 1);
		expr.end();
	}
	// second constraints set
	IloNumVar z(env, 0, IloInfinity, "ObjZ");
	for (const auto q : allPositions)
	{
		IloExpr expr(env);
		for (const auto it : t_allItems)
		{
			// calculate W(j, q)
			for (const auto& it2 : t_mapPosWidth.find(it->idx)->second)
			{
				if (it2 <= q && it2 >= q - it->width + 1)
				{
					auto iter = allVars.find(StripPacking::getVarName(it->idx, it2));
					assert(iter != allVars.end());
					expr += iter->second*it->height;
				}
			}
		}
		model.add(expr <= z);
		expr.end();
	}
	model.add(IloMinimize(env, z));
	IloCplex cplex(env);
	cplex.extract(model);
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	//cplex.exportModel("lowerBound5.lp");
	
	cplex.setParam(IloCplex::Param::Preprocessing::RepeatPresolve, 3);
	cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 3);
	cplex.setParam(IloCplex::Param::MIP::Strategy::Probe, 3);
	cplex.setParam(IloCplex::Param::Preprocessing::Symmetry, 5);
	cplex.solve();
	double result = cplex.getObjValue();
	env.end();
	return result;
}



int StripPacking::subSetSum(const std::vector<int>& t_v, const int t_limit)
/*
Args:
	Given a vector of integers and a limit. Find the subset of the integers of which the sum is the largest but not exceeding the limit
	This is achieved by a simple dynamic programming algorithm
Returns:
	The best sum.
*/
{
	int** values = new int*[t_v.size() + 1];
	for (int i = 0; i < t_v.size() + 1; ++i)
	{
		values[i] = new int[t_limit + 1];
		if (i == 0)
		{
			for (int j = 0; j < t_limit + 1; ++j) values[i][j] = 0;
		}
		values[i][0] = 0;
	}

	for (int i = 1; i < t_v.size() + 1; ++i)
	{
		for (int j = 1; j < t_limit + 1; ++j)
		{
			if (t_v[i - 1] > j) values[i][j] = values[i - 1][j];
			else values[i][j] = std::max(values[i - 1][j - t_v[i - 1]] + t_v[i - 1], values[i - 1][j]);
		}
	}
	int result = values[t_v.size()][t_limit];
	for (int i = 0; i < t_v.size(); ++i) delete values[i];
	delete values;
	return result;
}
