#pragma once
#include <vector>

constexpr double tolerance = 0.0001;

template<typename T>
int dynamicPrg4KnapSack(const std::vector<const T*>& items, int capacity)
{
	size_t itemSize = items.size() + 1;
	size_t weightSize = capacity + 1;
	int** valueMatrix = new int*[weightSize];
	for (size_t i = 0; i < weightSize; ++i)
	{
		valueMatrix[i] = new int[itemSize];
		for (size_t j = 0; j < itemSize; ++j)
		{
			valueMatrix[i][j] = 0;
		}
	}
	for (size_t j = 1; j < itemSize; ++j)
	{
		for (size_t i = 1; i < weightSize; ++i)
		{
			if (items[j - 1]->weight > i)
			{
				valueMatrix[i][j] = valueMatrix[i][j - 1];
			}
			else
			{
				valueMatrix[i][j] = std::max((valueMatrix[(i - items[j - 1]->weight)][j - 1] + items[j - 1]->value), (valueMatrix[i][j - 1]));
			}
		}
	}
	int result = valueMatrix[weightSize - 1][itemSize - 1];
	for (size_t i = 0; i < weightSize; ++i)
	{
		delete[] valueMatrix[i];
	}
	delete[] valueMatrix;
	return result;
}


/*
return the last row of the valueMatrix
*/
template<typename T>
std::vector<int> dynamicPrg4KnapSack(const std::vector<const T*>& items, int capacity, int dummy)
{
	std::vector<int> result(items.size(), 0);
	size_t itemSize = items.size() + 1;
	size_t weightSize = capacity + 1;
	int** valueMatrix = new int*[weightSize];
	for (size_t i = 0; i < weightSize; ++i)
	{
		valueMatrix[i] = new int[itemSize];
		for (size_t j = 0; j < itemSize; ++j)
		{
			valueMatrix[i][j] = 0;
		}
	}
	for (size_t j = 1; j < itemSize; ++j)
	{
		for (size_t i = 1; i < weightSize; ++i)
		{
			if (items[j - 1]->weight > i)
			{
				valueMatrix[i][j] = valueMatrix[i][j - 1];
			}
			else
			{
				valueMatrix[i][j] = std::max((valueMatrix[(i - items[j - 1]->weight)][j - 1] + items[j - 1]->value), (valueMatrix[i][j - 1]));
			}
		}
	}
	for (size_t i = 1; i <= itemSize - 1; ++i)
		result[i-1] = valueMatrix[weightSize - 1][i];
	for (size_t i = 0; i < weightSize; ++i)
	{
		delete[] valueMatrix[i];
	}
	delete[] valueMatrix;
	return result;
}

double dynamicPrg4KnapSack(const std::vector<double>& t_values, const std::vector<int>& t_weights, int capacity, std::vector<int>& t_selected);
