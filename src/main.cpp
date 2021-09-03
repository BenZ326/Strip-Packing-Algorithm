#include "datareader.h"
#include <experimental/filesystem>
#include <string>
#include "spp.h"
#include <iostream>
#include <map>
#include "BLEU.h"
#include "heuristic.h"
int main()
{
	const std::string instancesFolder = "./2sp/";
	for (int i = 0; i < 1; ++i) {
	for (const auto& entry : std::experimental::filesystem::directory_iterator(instancesFolder))
	{
		std::string filePath = entry.path().relative_path().string();
		std::cout << filePath;
		std::vector<const StripPacking::item*> allItems;
		StripPacking::Heuristic hrs;
		int W = readData(filePath, allItems);
		std::vector<const StripPacking::item*> copyItems(allItems.begin(), allItems.end());
		int totalArea = 0;
		StripPacking::BLEU alg(allItems,W, 20,1000);
		auto status = alg.evaluate();
		std::cout << "The status is " << status<<"\n";
		for (auto it = allItems.begin(); it != allItems.end(); ++it)
			delete (*it);
	}
	}
	system("pause");
}