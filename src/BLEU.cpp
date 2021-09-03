#include "BLEU.h"
#include <algorithm>
#include "spp.h"
#include "knapsack.h"
#include <iostream>
#include <stack>

#include <chrono>
double StripPacking::BLEU::tolerance = 0.0001;
int StripPacking::BLEU::bigNumber = 999999;
int StripPacking::BLEU::BBMaxExplNodesPerPack = 10000000;
int StripPacking::BLEU::BBMaxExplNodesNonPerPack = 80000;
int StripPacking::BLEU::interestingStatics = 0;
int StripPacking::BLEU::ycheckExplNode = 10000000;
bool StripPacking::BLEU::nodeLimitFlag = false;
StripPacking::algorithmStatus StripPacking::BLEU::algStatus = StripPacking::algorithmStatus::exact;
/*
only invoke when it's in evaluatedMode
*/
StripPacking::BLEU::BLEU(const std::vector<const item*>& t_items, const int t_W, const int t_TrialHeight,
	const int t_timeLimit)
	:_allItems(t_items),_W(t_W), _trialHeight(t_TrialHeight), _evaluatedMode(true), _timeLimit(t_timeLimit)
{
	// sort the items by the nonincreasing of width and breaking ties by nonincreasing height
	std::sort(_allItems.begin(), _allItems.end(), compareItemByWidth);
	this->reassignItemsIdx();
	this->preprocessing();
	this->bounds();
}

StripPacking::BLEU::BLEU(const std::vector<const item*>& t_items, const int t_W, const int t_timeLimit)
	:_allItems(t_items), _W(t_W), _evaluatedMode(false), _timeLimit(t_timeLimit)
{
	// sort the items by the nonincreasing of width and breaking ties by nonincreasing height
	std::sort(_allItems.begin(), _allItems.end(), compareItemByWidth);
	this->reassignItemsIdx();
	this->preprocessing();
	this->bounds();
}

/*
reassign idx for items
*/
void StripPacking::BLEU::reassignItemsIdx()
{
	int idx = 0;
	for (auto it = _allItems.begin(); it != _allItems.end(); ++it)
	{
		const_cast<item*>(*it)->idx = idx++;
	}
}



const StripPacking::solutionStatus StripPacking::BLEU::evaluate()
{
	if (_processedItems.empty()) return solutionStatus::feasible;
	if (_bestLowerBound > _trialHeight) return solutionStatus::infeasible;
	int binWidth = _processedW;
	int binHeight = _trialHeight - _processedH;
	int increment = 0;
	std::vector<item*> tmpItems;
	for (const auto& it : _processedItems)
	{
		item* tmp = new item(it->idx, it->width, it->height, it->idxHelper);
		tmpItems.push_back(tmp);
	}
	this->preprocessItemHeight(tmpItems, binHeight, binWidth);
	if (binWidth < 0) return StripPacking::solutionStatus::infeasible;
	// make items constant
	std::vector<const item*> Items;
	for (const auto& it : tmpItems) Items.push_back(std::move(it));
	// end preprocess
	auto status = this->branchAndBound(Items, binWidth, binHeight);
	BLEU::algStatus = (BLEU::nodeLimitFlag && status == solutionStatus::infeasible) ? algorithmStatus::approximate : BLEU::algStatus;
	this->releaseTmpItems(Items);
	return status;
}


const StripPacking::solutionStatus StripPacking::BLEU::solvePCC()		
// solve the parallel machine scheduling with contiguity constraints
{
	if (_processedItems.empty()) return solutionStatus::feasible;
	if (_bestLowerBound > _trialHeight) return solutionStatus::infeasible;
	int binWidth = _processedW;
	int binHeight = _trialHeight - _processedH;
	int increment = 0;
	std::vector<item*> tmpItems;
	for (const auto& it : _processedItems)
	{
		item* tmp = new item(it->idx, it->width, it->height, it->idxHelper);
		tmpItems.push_back(tmp);
	}
	this->preprocessItemHeight(tmpItems, binHeight, binWidth);
	if (binWidth < 0) return StripPacking::solutionStatus::infeasible;
	// make items constant
	std::vector<const item*> Items;
	for (const auto& it : tmpItems) Items.push_back(std::move(it));
	// end preprocess
	auto status = this->branchAndBoundYRelax(Items, binWidth, binHeight);
	this->releaseTmpItems(Items);
	return status;

}


/*
the y-check algorithm ---------------------------------------------------------------------------------------start
*/

/*
Arguements:
t_processedW: the width of the strip
t_TrialHeight: the tempting height
itemPositions: the left-bottom coordinates of all the items
t_processedItems: all the items to be packed in the strip
Return: if it is a feasible solution for the SPP then return true, else false
*/
bool StripPacking::BLEU::yCheckAlgorithm(const int t_processedW, const int t_TrialHeight, const std::vector<coordinate>& itemPositions,
	const std::vector<const item*> t_processedItems) const

{
	std::vector<coordinate> Cords4yCheck = itemPositions;
	int binWidth = t_processedW;
	auto Items = this->preprocess4yCheck(binWidth, t_processedItems, Cords4yCheck, t_TrialHeight);
	bool result = (this->yCheckEnumerationTree(Items, Cords4yCheck, t_TrialHeight, binWidth) == solutionStatus::feasible);
	this->releaseTmpItems(Items);
	//bool result = (this->yCheckEnumerationTree(t_processedItems, itemPositions, t_TrialHeight
	//	, t_processedW) == solutionStatus::feasible);
	return result;
}

/*

t_InterestItems: the items that are going through the enumeration tree
t_xCords: the corresonding x coordinates of the items, with the same ordering as that of t_interestItems
t_Height: the height of a bin to pack items in t_InterestItems
t_Width: the width of a bin to pack items in t_InterestItems
*/
StripPacking::solutionStatus StripPacking::BLEU::yCheckEnumerationTree(const std::vector<const item*>& t_InterestItems, const std::vector<coordinate>& t_Cords,
	const int t_Height, const int t_Width) const
{
	if (t_Width == 0) return solutionStatus::infeasible;			// it could happen when the very item ends at the current last column, then there will be no space for any merge
	std::unique_ptr<BBNode> root(new BBNode(t_InterestItems, t_Cords, t_Width, t_Height));
	std::stack<std::unique_ptr<BBNode>> yEnTree;
	yEnTree.push(std::move(root));
	int exploreNodes = 0;
	while (!yEnTree.empty())
	{
		const auto currentNode = std::move(yEnTree.top());
		yEnTree.pop();
		if (currentNode->remainingItems.empty())
		{
			return solutionStatus::feasible;
		}
		if (this->yCheckBounding(currentNode))	continue;
		exploreNodes++;
		this->yCheckMakeBranch(currentNode, yEnTree);
		if (exploreNodes > StripPacking::BLEU::ycheckExplNode)
		{
			StripPacking::BLEU::nodeLimitFlag = true;
			return solutionStatus::pending;
		}
	}
	return solutionStatus::infeasible;
}

bool StripPacking::BLEU::yCheckBounding(const std::unique_ptr<BBNode>& t_currentNode) const
{
	// fathoming criteria 1
	for (size_t i = 0; i < t_currentNode->columnsOccupiedHeight.size(); ++i)
	{
		int sumHeight = t_currentNode->columnsOccupiedHeight[i];
		for (auto it : t_currentNode->remainingItems)
		{
			if (t_currentNode->itemPositions[it->idxHelper].x <= i &&
				t_currentNode->itemPositions[it->idxHelper].x + it->width-1 >= i)
			{
				sumHeight += it->height;
			}
			if (sumHeight > t_currentNode->trialHeight)
				return true;
		}
	}

	if (!t_currentNode->packedItems.empty())
	{
		const item* tmpItem = t_currentNode->packedItems.back();
		// fathoming criteria 3
		for (const auto& it : t_currentNode->remainingItems)
		{
			if (tmpItem->height == it->height && tmpItem->width == it->width &&
				t_currentNode->itemPositions[tmpItem->idxHelper].x == t_currentNode->itemPositions[it->idxHelper].x&&
				it->idx < tmpItem->idx)
				return true;
		}
	}
	
	if (!t_currentNode->packedItems.empty()) {
		// fathoming criteria 4
		const item* itemJ = t_currentNode->packedItems.back();
		auto p_js = t_currentNode->itemPositions[itemJ->idxHelper].x;
		for (size_t i = 0; i < t_currentNode->packedItems.size() - 1; ++i)
		{
			const item* itemK = t_currentNode->packedItems[i];
			auto xCord = t_currentNode->itemPositions[itemK->idxHelper].x;
			if (p_js == xCord && itemJ->idx < itemK->idx && itemJ->width == itemK->width)
			{
				if (t_currentNode->itemPositions[itemK->idxHelper].y == t_currentNode->columnsOccupiedHeight[p_js] - itemJ->height- itemK->height)
				{
					return true;
				}
			}
		}
	}
	return false;
}

void StripPacking::BLEU::yCheckMakeBranch(const std::unique_ptr<BBNode>& t_currentNode,
	std::stack<std::unique_ptr<BBNode>>& t_yEntree) const
{
	std::list<std::unique_ptr<BBNode>> children;
	std::vector<const item*> copyRemainingItems = t_currentNode->remainingItems;
	std::sort(copyRemainingItems.begin(), copyRemainingItems.end(), compareItemByxCords(t_currentNode->itemPositions));
	auto niche = this->getNiche(t_currentNode->columnsOccupiedHeight);
	// h(l)
	int l_niche = niche[0] == 0 ? t_currentNode->trialHeight : t_currentNode->columnsOccupiedHeight[niche[0] - 1];
	// h(r)
	int r_niche = niche[1] == t_currentNode->columnsOccupiedHeight.size() - 1 ? t_currentNode->trialHeight : t_currentNode->columnsOccupiedHeight[niche[1] + 1];
	bool emptyItem = true;
	for (size_t i = 0; i < copyRemainingItems.size(); ++i)
	{
		if (t_currentNode->itemPositions[copyRemainingItems[i]->idxHelper].x >= niche[0] &&
			t_currentNode->itemPositions[copyRemainingItems[i]->idxHelper].x + copyRemainingItems[i]->width - 1 <= niche[1])
		{
			bool fathom = false;
			const int& p_js = t_currentNode->itemPositions[copyRemainingItems[i]->idxHelper].x;
			 //fathoming criteria 5
			if (p_js > niche[0])
			{
				for (size_t k= 0; k < copyRemainingItems.size(); ++k)
				{
					if (t_currentNode->itemPositions[copyRemainingItems[k]->idxHelper].x >= niche[0] &&
						t_currentNode->itemPositions[copyRemainingItems[k]->idxHelper].x + copyRemainingItems[k]->width - 1 <= niche[1])
					{
						if (i == k)
							continue;
						const int& p_ks = t_currentNode->itemPositions[copyRemainingItems[k]->idxHelper].x;
						if (p_ks + copyRemainingItems[k]->width -1 <= p_js)
						{
							if (copyRemainingItems[k]->height <= std::min(l_niche - t_currentNode->columnsOccupiedHeight[niche[0]], copyRemainingItems[i]->height))
							{
								fathom = true;
								break;
							}
						}
					
					}
				}
			
			}
			if (fathom) continue;
			std::unique_ptr<BBNode> child(new BBNode(*t_currentNode, copyRemainingItems[i]));
			child->itemPositions[copyRemainingItems[i]->idxHelper].y = 
				child->columnsOccupiedHeight[child->itemPositions[copyRemainingItems[i]->idxHelper].x];
			for (size_t j = child->itemPositions[copyRemainingItems[i]->idxHelper].x; 
				j < child->itemPositions[copyRemainingItems[i]->idxHelper].x + copyRemainingItems[i]->width; ++j)
			{
				child->columnsOccupiedHeight[j] += copyRemainingItems[i]->height;
			}
			// fathoming criteria 2
			if (emptyItem && child->columnsOccupiedHeight[p_js] + copyRemainingItems[i]->height<= std::min(l_niche, r_niche))
				emptyItem = false;
			children.push_back(std::move(child));
		}
	}

	if (emptyItem)
	{
		std::unique_ptr<BBNode> child(new BBNode(*t_currentNode));
		for (size_t j = niche[0]; j <= niche[1]; ++j)
		{
			child->columnsOccupiedHeight[j] = std::min(l_niche, r_niche);
		}
		children.push_back(std::move(child));
	}
	for (auto it = children.begin(); it != children.end(); ++it)
		t_yEntree.push(std::move(*it));
}
/*
Get the Niche 
t_ColumnsHeight: it is the height of each column 
the return value is a two dimensional array, one element of the array is the start column of the niche
while the second is the end column of the niche
*/
std::vector<int> StripPacking::BLEU::getNiche(const std::vector<int>& t_ColumnHeights) const
{
	std::vector<int>result (2, 0);
	auto ptr = std::min_element(t_ColumnHeights.begin(), t_ColumnHeights.end());
	int startCol =  ptr - t_ColumnHeights.begin();
	int endCol = t_ColumnHeights.size()-1;
	for (size_t i = startCol+1; i < t_ColumnHeights.size(); i++)
	{
		if (t_ColumnHeights[i] > (*ptr))
		{
			endCol = i - 1;
			break;
		}
	}
	result[0] = startCol;
	result[1] = endCol;
	return result;
}

/*
the y-check algorithm ---------------------------------------------------------------------------------------end
*/


//The branch and bound algorithm -----------------------------------------------------------------------------start
/*

*/
StripPacking::BLEU::BBNode::BBNode(const std::vector<const item*>& t_remainingItems, const int t_Width, const int t_TrialHeight):
remainingItems(t_remainingItems), trialHeight(t_TrialHeight)
{
	leftMostIdx = 0;
	columnsOccupiedHeight = std::vector<int>(t_Width,0);			// the order is consistent to the left most column -> the right most column
	maxiItemIdxColumns = std::vector<int>(t_Width,-1);				// the order is consistent to the left most column -> the right most column
	coordinate dummy(-1, -1);
	itemPositions = std::vector<coordinate>(t_remainingItems.size(), dummy);
}

StripPacking::BLEU::BBNode::BBNode(const BBNode& t_BBNode):trialHeight(t_BBNode.trialHeight)
{
	leftMostIdx = t_BBNode.leftMostIdx;
	columnsOccupiedHeight = t_BBNode.columnsOccupiedHeight;			// [10,5,3,2] means 10 units of height in the 1st column is occupied and 5 units for the 2nd column...
	remainingItems = t_BBNode.remainingItems;			// make heap
	maxiItemIdxColumns = t_BBNode.maxiItemIdxColumns;				// [3,5,1,7] means among items placed in 1st column, 3 is the largest index of...
	itemPositions = t_BBNode.itemPositions;
	packedItems = t_BBNode.packedItems;
}

// build a BBNode for the y check algorithm
StripPacking::BLEU::BBNode::BBNode(const std::vector<const item*>& t_remainingItems, const std::vector<coordinate>& t_Cords,
	const int t_Width, const int t_TrialHeight):remainingItems(t_remainingItems), trialHeight(t_TrialHeight)
{
	leftMostIdx = 0;
	columnsOccupiedHeight = std::vector<int>(t_Width, 0);			// the order is consistent to the left most column -> the right most column
	maxiItemIdxColumns = std::vector<int>(t_Width, 0);				// the order is consistent to the left most column -> the right most column
	coordinate dummy(-1, -1);
	itemPositions = t_Cords;
}

StripPacking::BLEU::BBNode::BBNode(const BBNode& t_BBNode, const item* t_placedItem)			// invoked when making a branch, the t_placedItem is the packed item in this time of branching
:trialHeight(t_BBNode.trialHeight)
{
	packedItems = t_BBNode.packedItems;
	packedItems.push_back(t_placedItem);
	for (const auto& it : t_BBNode.remainingItems)
	{
		if (it->idx == t_placedItem->idx)
			continue;
		else this->remainingItems.push_back(it);
	}
	leftMostIdx = t_BBNode.leftMostIdx;
	columnsOccupiedHeight = t_BBNode.columnsOccupiedHeight;			// [10,5,3,2] means 10 units of height in the 1st column is occupied and 5 units for the 2nd column...
	maxiItemIdxColumns = t_BBNode.maxiItemIdxColumns;				// [3,5,1,7] means among items placed in 1st column, 3 is the largest index of...
	itemPositions = t_BBNode.itemPositions;
}

const bool StripPacking::BLEU::bounding(const std::unique_ptr<BBNode>& t_currentNode) const
{
	//fathoming criteria 1
	if (!t_currentNode->packedItems.empty())
	{
		const item* tmpItem = t_currentNode->packedItems.back();
		for (const auto & it : t_currentNode->remainingItems)
		{
			if (it->height == tmpItem->height && it->width == tmpItem->width && it->idx < tmpItem->idx)
				return true;
		}
	}
	// fathoming criteria 2 is merged in the makeBranch function

	// standard continuous bounding which is described in the section 5.2 "branch and bound for the spp(L)"
	// fathoming criteria 3
	int remainingArea = 0;
	for(const auto& it : t_currentNode->remainingItems) 	remainingArea += it->height*it->width;
	int spaceArea = 0;
	for (size_t i = 0; i < t_currentNode->columnsOccupiedHeight.size(); ++i)
		spaceArea += (t_currentNode->trialHeight - t_currentNode->columnsOccupiedHeight[i]);
	if (remainingArea > spaceArea)
		return true;
	// fathoming criteria 4
	 //dynamic cuts:
	if (this->dynamicCuts(t_currentNode))
	{
		BLEU::interestingStatics++;
		return true;
	}
	return false;
	
}

void  StripPacking::BLEU::makeBranch(const std::unique_ptr<BBNode>& t_currentNode,
 std::stack<std::unique_ptr<BBNode>>& t_dfstree) const
{
	std::set<whPair> packedPairs;
	std::list<std::unique_ptr<BBNode>> children;
	std::list<std::unique_ptr<BBNode>> emptyChild;
	// selects the left-most column
	int selectedColumn = t_currentNode->leftMostIdx;   // start from 0, so it ranges from 0, 1,2,3,....., _processedW-1
	// pack nothing
	if (t_currentNode->columnsOccupiedHeight[selectedColumn] > 0)
	{
		std::unique_ptr<BBNode> child(new BBNode(*t_currentNode.get()));
		child->columnsOccupiedHeight[selectedColumn] = child->trialHeight;
		child->leftMostIdx++;
		emptyChild.push_back(std::move(child));
	}
	// pack j on the column
	std::vector<const item*> copyRemainingItems = t_currentNode->remainingItems;
	std::sort(copyRemainingItems.begin(), copyRemainingItems.end(),compareItemByIdx);
	for(int i = 0;  i< copyRemainingItems.size();i++)
	{
		auto chosenItem = copyRemainingItems[i];
		// check width of the item 
		if (chosenItem->width + selectedColumn > t_currentNode->columnsOccupiedHeight.size())
			continue;
		// check height of the item
		if (chosenItem->height + t_currentNode->columnsOccupiedHeight[selectedColumn] > t_currentNode->trialHeight)
			continue;
		if (t_currentNode->maxiItemIdxColumns[selectedColumn] > chosenItem->idx)
			continue;
		int minHeight = 99999;
		for (size_t k = 0; k < copyRemainingItems.size(); ++k)
		{
			if (k == i) continue;
			if (minHeight > copyRemainingItems[k]->height)
				minHeight = copyRemainingItems[k]->height;
		}
		// make a child by pack the item on the column, update the coordinate 
		std::unique_ptr<BBNode> child(new BBNode(*t_currentNode.get(), copyRemainingItems[i]));
		child->itemPositions[chosenItem->idxHelper].x = selectedColumn;
		child->itemPositions[chosenItem->idxHelper].y = t_currentNode->columnsOccupiedHeight[selectedColumn];
		child->maxiItemIdxColumns[selectedColumn] = chosenItem->idx;
		//if there exists an item can be added on the column
		for (int offset = 0; offset <= chosenItem->width - 1; ++offset)
		{
			if (chosenItem->height + t_currentNode->columnsOccupiedHeight[selectedColumn + offset] + minHeight
				<= t_currentNode->trialHeight)
			{
				child->columnsOccupiedHeight[selectedColumn + offset] += chosenItem->height;
			}
			else
			{
				child->columnsOccupiedHeight[selectedColumn + offset] = t_currentNode->trialHeight;
				child->leftMostIdx = selectedColumn + offset +1;
			}
		}
		children.push_back(std::move(child));
	}
	// load the empty child
	if (!emptyChild.empty())
		t_dfstree.push(std::move(emptyChild.back()));
	for (auto it = children.rbegin(); it != children.rend(); ++it)
		t_dfstree.push(std::move((*it)));
}

const StripPacking::solutionStatus StripPacking::BLEU::branchAndBound(const std::vector<const item*>& t_Items, const int t_binWidth,
	const int t_binHeight)
{
	if (t_Items.empty()) return StripPacking::solutionStatus::feasible;
	StripPacking::BLEU::nodeLimitFlag = false;
	BLEU::algStatus = algorithmStatus::exact;
	int tmpW = t_binWidth;
	int tmpH = t_binHeight;
	BLEU::interestingStatics = 0;
	int maxExpNodes;
	std::unique_ptr<BBNode>	root(new BBNode(t_Items, tmpW, tmpH));
	double totalArea = 0.0;
	for (const auto& it : t_Items) totalArea += it->width*it->height;
	if (abs(tmpH-(totalArea / tmpW)) <BLEU::tolerance) maxExpNodes = BLEU::BBMaxExplNodesPerPack;
	else maxExpNodes = BLEU::BBMaxExplNodesNonPerPack;
	std::stack<std::unique_ptr<BBNode>> dfsTree;
	dfsTree.push(std::move(root));
	int numberExploredNodes = 0;
	while (!dfsTree.empty() && numberExploredNodes < maxExpNodes)
	{
		const auto currentNode = std::move(dfsTree.top());
		dfsTree.pop();
		// if it's a feasible solution then invoke the y-check algorithm
		if (currentNode->remainingItems.empty())
		{
			if (this->yCheckAlgorithm(tmpW, tmpH, currentNode->itemPositions, t_Items))
			{
				return solutionStatus::feasible;
			}
			else continue;							// the node can not be transformed to a feasible solution for the SPP
		}
		else
		{
			// bounding the current Node
			if (this->bounding(currentNode))
				continue;
			// make branch
			numberExploredNodes++;
			this->makeBranch(currentNode, dfsTree);
		}
	}
	if (numberExploredNodes >= maxExpNodes)
	{
		return solutionStatus::pending;
	}
	else return solutionStatus::infeasible;
}


const StripPacking::solutionStatus StripPacking::BLEU::branchAndBoundYRelax(const std::vector<const item*>& t_Items, const int t_binWidth,
	const int t_binHeight)
{
	if (t_Items.empty()) return StripPacking::solutionStatus::feasible;
	StripPacking::BLEU::nodeLimitFlag = false;
	int tmpW = t_binWidth;
	int tmpH = t_binHeight;
	BLEU::interestingStatics = 0;
	int maxExpNodes;
	std::unique_ptr<BBNode>	root(new BBNode(t_Items, tmpW, tmpH));
	double totalArea = 0.0;
	for (const auto& it : t_Items) totalArea += it->width*it->height;
	if (abs(tmpH - (totalArea / tmpW)) < BLEU::tolerance) maxExpNodes = BLEU::BBMaxExplNodesPerPack;
	else maxExpNodes = BLEU::BBMaxExplNodesNonPerPack;
	std::stack<std::unique_ptr<BBNode>> dfsTree;
	dfsTree.push(std::move(root));
	int numberExploredNodes = 0;
	while (!dfsTree.empty() && numberExploredNodes < maxExpNodes)
	{
		const auto currentNode = std::move(dfsTree.top());
		dfsTree.pop();
		// if it's a feasible solution return 
		if (currentNode->remainingItems.empty())
		{
			return solutionStatus::feasible;
		}
		else
		{
			// bounding the current Node
			if (this->bounding(currentNode))
				continue;
			// make branch
			numberExploredNodes++;
			this->makeBranch(currentNode, dfsTree);
		}
	}
	if (numberExploredNodes >= maxExpNodes)
	{
		return solutionStatus::pending;
	}
	else return solutionStatus::infeasible;
}

const bool StripPacking::BLEU::dynamicCuts(const std::unique_ptr<BBNode>& t_currentNode) const
{
	std::list<coordinate> leftCorners;
	int prev = t_currentNode->trialHeight;
	for (size_t i=0; i<t_currentNode->columnsOccupiedHeight.size(); ++i)
	{
		if (prev > t_currentNode->columnsOccupiedHeight[i])
		{
			coordinate cords(i, t_currentNode->columnsOccupiedHeight[i]);
			leftCorners.push_back(cords);
		}
		prev = t_currentNode->columnsOccupiedHeight[i];
	}
	bool result = false;
	size_t tmpSize = leftCorners.size();
	/*
	s[i][0]: the width of the i^th stage,
	s[i][1]: the height of the i^th stage
	g[i][0]: the x coordinate of the i^th left corner
	g[i][1]: the y coordinate of the i^th left corner
	l[i][0]: the realizable maximal width
	l[i][1]: the realizable maximal height
	*/
	int** s = new int*[tmpSize];
	int** g = new int*[tmpSize];
	int** l = new int*[tmpSize];
	for (size_t i = 0; i < tmpSize; ++i)
	{
		s[i] = new int[2];
		g[i] = new int[2];
		l[i] = new int[2];
	}
	int i = 0;
	int tmpXPrev = -1;			// store the x coordinate of the previous left corner
	int tmpYPrev = -1;			// store the y coordinate of the previous left corner
	/*
	At each left corner, the length can be updated for the onging left corner,
	while the width can only be updated by the upcoming left corner
	*/
	for (std::list<coordinate>::const_iterator iter = leftCorners.begin();
		iter != leftCorners.end(); ++iter)
	{
		g[i][0] = t_currentNode->maxiItemIdxColumns.size() - (*iter).x;
		g[i][1] = t_currentNode->trialHeight - (*iter).y;
		if (tmpXPrev == -1)
		{
			tmpXPrev = (*iter).x;
			tmpYPrev = (*iter).y;
			s[i][1] = t_currentNode->trialHeight - (*iter).y;
		}
		else
		{
			s[i - 1][0] = (*iter).x - tmpXPrev;
			s[i][1] = tmpYPrev - (*iter).y;
			tmpXPrev = (*iter).x;
			tmpYPrev = (*iter).y;
		}
		if (i == tmpSize - 1)
			s[i][0] = t_currentNode->maxiItemIdxColumns.size() - (*iter).x;
		++i;
	}

	/*
	compute the remaining available area (in the paper, it is denoted as V_pi)
	*/
	int vPi = 0;
	for (size_t i = 0; i < tmpSize; ++i)
		vPi += g[i][1] * s[i][0];

	/*
	preapre the onedimensional items for running the DP
	*/
	std::vector<const oneDimensionItem*> oneDim4HeightVec;
	std::vector<const oneDimensionItem*> oneDim4WidthVec;
	int accumulatedArea = 0;			//accumulated Area of the items considered so far
	for (std::vector<const item*>::const_iterator iter = t_currentNode->remainingItems.begin();
		iter != t_currentNode->remainingItems.end(); ++iter)
	{
		
		const oneDimensionItem* oneDimItemH = new oneDimensionItem((*iter)->height, (*iter)->height);
		const oneDimensionItem* oneDimItemW = new oneDimensionItem((*iter)->width, (*iter)->width);
		oneDim4HeightVec.push_back(oneDimItemH);
		oneDim4WidthVec.push_back(oneDimItemW);
	}
	std::vector<std::vector<int>> l0ValueOverItems;
	std::vector<std::vector<int>> l1ValueOverItems;
	for (size_t i = 0; i < tmpSize; ++i)
	{
		l0ValueOverItems.push_back(dynamicPrg4KnapSack<oneDimensionItem>(oneDim4WidthVec, g[i][0],0));
		l1ValueOverItems.push_back(dynamicPrg4KnapSack<oneDimensionItem>(oneDim4HeightVec, g[i][1],0));
	}
	/*
		calculate the values for matrix l
		*/
	for (size_t j = 0; j<t_currentNode->remainingItems.size();++j)
	{
		accumulatedArea += t_currentNode->remainingItems [j]->height*t_currentNode->remainingItems[j]->width;
		int A = 0;
		int sum_B = 0;
		std::vector<int> B;
		B.clear();
		for (size_t i = 0; i < tmpSize; ++i)
		{
			l[i][0] = l0ValueOverItems[i][j];
			l[i][1] = l1ValueOverItems[i][j];
		}
		for (size_t i = 0; i < tmpSize; ++i)
		{
			A += s[i][0] * (g[i][1] - l[i][1]) + s[i][1] * (g[i][0] - l[i][0]) - (g[i][1] - l[i][1])*(g[i][0] - l[i][0]);
			if (i < tmpSize - 1)
			{
				A += (g[i][1] - l[i][1])*(g[i + 1][0] - l[i + 1][0]);
			}
			if (t_currentNode->remainingItems.size() < tmpSize)
			{
				if (i == 0)
					B.push_back((g[1][0] - l[1][0] + s[0][0] - g[0][0] + l[0][0])*(s[0][1] - g[0][1] + l[0][1]));
				if (i == tmpSize - 1)
					B.push_back((s[i][0] - g[i][0] + l[i][0])*(g[i - 1][1] - l[i - 1][1] + s[i][1] - g[i][1] + l[i][1]));
				if (i > 0 && i < tmpSize - 1)
					B.push_back((g[i + 1][0] - l[i + 1][0] + s[i][0] - g[i][0] + l[i][0])*(g[i - 1][1] - l[i - 1][1] + s[i][1] - g[i][1] + l[i][1]));
			}
		}
		if (!B.empty())
		{
			std::sort(B.begin(), B.end());
			for (size_t k = 0; k < leftCorners.size() - t_currentNode->remainingItems.size(); ++k)
			{
				sum_B += B[k];
			}
		}
		if (vPi - (A + sum_B) < accumulatedArea)
		{
			result = true;
			break;
		}
	}
		
	/*
	Release the heap
	*/
	{
		for (size_t i = 0; i < tmpSize; ++i)
		{
			delete[] s[i];
			delete[] g[i];
			delete[] l[i];
		}
		delete[] s;
		delete[] g;
		delete[] l;
		for (std::vector<const oneDimensionItem*>::iterator iter = oneDim4HeightVec.begin(); iter != oneDim4HeightVec.end();
			++iter)
			delete (*iter);
		for (std::vector<const oneDimensionItem*>::iterator iter = oneDim4WidthVec.begin(); iter != oneDim4WidthVec.end();
			++iter)
			delete (*iter);
	}
	return result;
}


//The branch and bound algorithm -----------------------------------------------------------------------------end


/*
Cut all items into pieces along the direction of height
*/
void StripPacking::BLEU::cutItemsAlongHeight()
{
	int idx = 0;
	for (const auto& it : _allItems)
	{
		for (size_t i = 1; i <= it->height; ++i)
		{
			const itemPieceWidth* tmp = new itemPieceWidth(idx++, it->width);
			_allItemPiecesWidths.push_back(tmp);
		}
	}
}

/*
5.1 the three preprocessing procedures
*/
void StripPacking::BLEU::preprocessing()
{
	this->preprocessingFixItems();
	this->preprocessingReduceW();
	this->preprocessingModifyItemWidth();
}

/*
An exact algorithm for the two-dimensional strip packing problem 2010 by Marco Antonio Boschetti, Lorenza Montaletti
Section 2.2.1
find a subset of items containing all items that cannot be placed side by side with any other items
*/
void StripPacking::BLEU::preprocessingFixItems()
{
	// find the item with the minimal width
	int minWidth = _allItems.back()->width;
	int sumHeight = 0;
	int idx = 0;
	for (const auto& it : _allItems)
	{
		if (it->width + minWidth > _W)
			sumHeight += it->height;
		else
		{
			const_cast<item*>(it)->idxHelper = idx++;				// idxHelper stores the idx in the _processedItems.
			_processedItems.push_back(it);
		}
	}
	_processedH = sumHeight;				// store the occupied height
}

/*
Calculate the maximal width that a subset of items can be packed side by side if the maximal width is strictly less than W,
than W = the maximal width
*/
void StripPacking::BLEU::preprocessingReduceW()
{
	for (const auto& it : _processedItems)
		_allWidths.push_back(it->width);
	_processedW = subSetSum(_allWidths,  _W);
}




/*
modify width for items according to 2.2.3 in the paper "An exact algorithm for the two-dimensional strip packing problem"
*/
void StripPacking::BLEU::preprocessingModifyItemWidth()
{
	for (std::vector<const item*>::iterator i_it = _processedItems.begin(); i_it != _processedItems.end(); ++i_it)
	{
		std::vector<int> integers;
		for (auto & j_it : _processedItems)
		{
			if ((*i_it)->idx == j_it->idx) continue;
			integers.push_back(j_it->width);
		}
		int maxWidth = subSetSum(integers, _processedW - (*i_it)->width) + (*i_it)->width;
		if (maxWidth < _processedW) 
			const_cast<item*>((*i_it))->width = (*i_it)->width + _processedW - maxWidth;
	}
}

/*
1) modify item height and preprocess t_items (delete a few items so the helper idx is also changed)
2) modify the bin width 
*/
void StripPacking::BLEU::preprocessItemHeight(std::vector<item*>& t_items, 
	const int t_binHeight, int& t_binWidth)
{
	for (auto& i_it : t_items)
	{
		std::vector<int> integers;
		for (auto & j_it : t_items)
		{
			if (i_it->idx == j_it->idx) continue;
			integers.push_back(j_it->height);
		}
		int maxHeight = subSetSum(integers, t_binHeight - i_it->height) + (i_it)->height;
		if (maxHeight < t_binHeight)
		{
			(i_it)->height = (i_it)->height + t_binHeight - maxHeight;
		}
	}
	int minHeight = BigNumber;
	for (const auto& it : t_items)
	{
		if (minHeight > it->height) 
			minHeight = it->height;
	}
	std::list<item*> tmpItems;
	for (const auto& it : t_items)
	{
		if (it->height + minHeight <= t_binHeight)
			tmpItems.push_back(it);
		else t_binWidth -= it->width;
	}
	t_items.clear();
	int i = 0;
	for (auto & it : tmpItems)
	{
		it->idxHelper = i;
		++i;
		t_items.push_back(it);
	}
}

/*
section 5.1 lower bounds plus upper bounds
*/
void StripPacking::BLEU::bounds()
{
	int lb1 = this->LowerBound1();
	int lb2 = this->LowerBound2();
	int lb4 = this->LowerBound4();
	_bestLowerBound = std::max({ lb1, lb2, lb4 });
	int lb3 = this->LowerBound3();
	int lb5 = this->LowerBound5();
	_bestLowerBound = std::max({ _bestLowerBound, lb3, lb5 });
}

const int StripPacking::BLEU::LowerBound1() const
{
	int sum = 0;
	for (const auto& it : _processedItems)
	{
		sum += it->height*it->width;
	}
	return int(std::ceil(sum / _W)) + _processedH;
}

/*
Section 3.2 in the paper "An exact algorithm for the two-dimensional strip packing problem" 
dual feasible functions 1, 2 ,3
For understanding the concept of dual feasible functions, refer to the paper : "New  classes of fast lower  bounds  for  bin  packingproblems"
*/
const int StripPacking::BLEU::LowerBound2() const
{
	//dual feasible function 1:
	int lowerBound = 0;
	for (size_t alpha = 1; alpha <= _processedW; ++alpha)
	{
		std::vector<double> allTransformedWidths;
		for (const auto& it : _processedItems) allTransformedWidths.push_back(this->DualFeasibleFunction1(alpha,it->width));
		double sum = 0.0;
		for (size_t i = 0; i < _processedItems.size(); i++)		sum += allTransformedWidths[i] * _processedItems[i]->height;
		if (lowerBound < std::ceil(sum / _processedW))
			lowerBound = std::ceil(sum / _processedW);
	}
	 //dual feasible function 2:
	for (size_t alpha = 1; alpha <= std::floor(_processedW/2); ++alpha)
	{
		std::vector<double> allTransformedWidths;
		for (const auto& it : _processedItems) 
			allTransformedWidths.push_back(this->DualFeasibleFunction2(alpha, it->width));
		double sum = 0.0;
		for (size_t i = 0; i < _processedItems.size(); i++)		sum += allTransformedWidths[i] * _processedItems[i]->height;
		auto candidate = std::ceil(sum / this->DualFeasibleFunction2(alpha, _processedW));
		if (lowerBound < candidate)
			lowerBound = candidate;
	}

	// dual feasible function 3:
	for (size_t alpha = 1; alpha <= std::floor(_processedW / 2); ++alpha)
	{
		std::vector<double> allTransformedWidths;
		for (const auto& it : _processedItems)
			allTransformedWidths.push_back(this->DualFeasibleFunction3(alpha, it->width));
		double sum = 0.0;
		for (size_t i = 0; i < _processedItems.size(); i++)		sum += allTransformedWidths[i] * _processedItems[i]->height;
		auto candidate = std::ceil(sum / this->DualFeasibleFunction3(alpha, _processedW));
		if (lowerBound < candidate)
			lowerBound = candidate;
	}
	return lowerBound + _processedH;
}

/*
The lower bound proposed in the paper:
Alvarez-Vales R, Parreno F, Tamarit JM. A branch and bound algorithm for the strip packing problem. OR spectrum. 
2009 Apr 1;31(2):431-59.
Section 4.2.6
*/
const int StripPacking::BLEU::LowerBound3() const
{
	// step 1:
	if (_processedItems.empty()) return _bestLowerBound;
	bool exitFlag = false;
	int result;
	for (int k = 0; !exitFlag; k++)
	{
		int RectangleW = _processedW;
		int RectangleH = _bestLowerBound + k;
		std::vector<item*> itemsArr;
		for (const auto & it : _processedItems)
		{
			item*  tmp = new item(*it);
			itemsArr.push_back(tmp);
		}
		while (true)
		{
			std::make_heap(itemsArr.begin(), itemsArr.end(), compareItemByWHDifference(RectangleH, RectangleW));
			std::pop_heap(itemsArr.begin(), itemsArr.end(), compareItemByWHDifference(RectangleH, RectangleW));
			item* selectedItem = itemsArr.back();
			itemsArr.pop_back();
			// check if the item can fit the rectangle 
			if (selectedItem->height > RectangleH && selectedItem->width > RectangleW) break;
			std::set<int> removedItems;
			// step 2 and step 3
			if (RectangleW - selectedItem->width <= RectangleH - selectedItem->height)
			{
				// pack the item at the bottom and update the rectangle and width of some items
				RectangleH -= selectedItem->height;
				for (size_t i = 0; i < itemsArr.size(); ++i)
				{
					if (itemsArr[i]->width <= RectangleW - selectedItem->width)
					{
						itemsArr[i]->height = std::max(0, itemsArr[i]->height - selectedItem->height);
						if (itemsArr[i]->height == 0)
							removedItems.insert(itemsArr[i]->idx);
					}
				}
			}
			else
			{
				// pack the item at the left
				RectangleW -= selectedItem->width;
				for (size_t i = 0; i < itemsArr.size(); ++i)
				{
					if (itemsArr[i]->height <= RectangleH - selectedItem->height)
					{
						itemsArr[i]->width = std::max(0, itemsArr[i]->width - selectedItem->width);
						if (itemsArr[i]->width == 0)
							removedItems.insert(itemsArr[i]->idx);
					}
				}
			}
			std::vector<item*> remainingItemsArr;
			for (size_t i = 0; i < itemsArr.size(); ++i)
			{
				if (removedItems.find(itemsArr[i]->idx) == removedItems.end()) 	remainingItemsArr.push_back(itemsArr[i]);
				else 	delete itemsArr[i];
			}
			if (remainingItemsArr.empty())   // means all the items are packed
			{
				exitFlag = true;
				result = _bestLowerBound + k;
				break;
			}
			itemsArr.clear();
			itemsArr.assign(remainingItemsArr.begin(), remainingItemsArr.end());
		}
	}
	return result;
}

/*
Section 5.1 
Solve a noncontiguous packing problem (NCBP)
*/
const int StripPacking::BLEU::LowerBound4()const
{
	std::vector<int> allWidths;
	int colIdx = 0;
	// initialize columns
	std::list<std::set<int>> columns;			// a collection of columns		a column is a set of item index
	for (size_t i=0; i< _processedItems.size(); ++i)
	{
		std::set<int> tmp;
		tmp.insert(_processedItems[i]->idx);
		columns.push_back(tmp);
		allWidths.push_back(_processedItems[i]->width);
	}
	// model the LP
	IloEnv env;
	IloModel NCBP(env);
	std::map<int, IloRange> allCstrs;
	// add constraints
	for (size_t i = 0; i < _processedItems.size(); ++i)
	{
		char nameCstr[64];
		sprintf_s(nameCstr, 64, "Cstr %d", _processedItems[i]->idx);
		IloExpr expr(env);
		IloRange cstr(env, _processedItems[i]->height, expr, IloInfinity, nameCstr);
		allCstrs.insert(std::pair<int, IloRange>(_processedItems[i]->idx, cstr));
		NCBP.add(cstr);
		expr.end();
	}
	// add the variables and the objective function
	IloExpr obj(env);
	for (const auto& i_it : columns)
	{
		IloNumColumn col(env);
		for (const auto& j_it : i_it)
		{
			auto ptr = allCstrs.find(j_it);
			if (ptr != allCstrs.end()) col += ptr->second(1);
		}
		char varName[64];
		sprintf_s(varName, 64, "Var %d", colIdx++);
		IloNumVar var(col);
		var.setName(varName);
		var.setBounds(0, +IloInfinity);
		NCBP.add(var);
		col.end();
		obj += var;
	}
	IloObjective realObj = IloMinimize(env, obj);
	NCBP.add(realObj);
	obj.end();
	try {
		IloCplex cplex(NCBP);
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		while (true)
		{
			cplex.solve();
			//cplex.exportModel("NCBP.lp");
			// solve the pricing problem
			std::vector<IloNum> dualValues;
			for (size_t i=0; i<_processedItems.size(); ++i)
				dualValues.push_back(cplex.getDual(allCstrs.find(_processedItems[i]->idx)->second));
			std::vector<int> selectedItems;
			double value = dynamicPrg4KnapSack(dualValues, allWidths, _processedW, selectedItems);
			if (1.0 - value < - BLEU::tolerance)		// negative reduced cost
			{
				// add a column
				std::set<int> newCol;
				for (const auto& it : selectedItems) newCol.insert(_processedItems[it]->idx);
				IloNumColumn col(env);
				for (const auto& it : newCol)
				{
					auto ptr = allCstrs.find(it);
					if (ptr != allCstrs.end()) 
						col += ptr->second(1);
				}
				char varName[64];
				sprintf_s(varName, 64, "Var %d", colIdx++);
				IloNumVar var(col);
				var.setName(varName);
				var.setBounds(0, +IloInfinity);
				NCBP.add(var);
				realObj.setLinearCoef(var, 1.0);
			}
			else
			{
				double objValue = cplex.getObjValue();
				env.end();
				double lowerBound;
				std::modf(objValue, &lowerBound);
				if (objValue - lowerBound > BLEU::tolerance)
					return lowerBound + 1 + _processedH;
				else
					return lowerBound + _processedH;	
			}
		}
		// add columns
		// return the solution
	}
	catch (IloException & e)
	{
		std::cout << "NCBP error\n";
		std::cout << e<< std::endl;
		return -1;
	}
	
}

/*
solve the root node of the parallel machine scheduling with contiguous constraints
*/
const int StripPacking::BLEU::LowerBound5()const
{
	int maxHeight = _bestLowerBound - _processedH;
	std::map<int, std::set<int>> mapPosWidth, mapPosHeight;
	for (size_t idx = 0; idx < _processedItems.size(); ++idx)
	{
		auto possiblePositionsWidth = computeFX(_processedW - _processedItems[idx]->width, idx,
			_processedItems, true);
		auto possiblePositionsHeight = computeFX(maxHeight - _processedItems[idx]->height, idx,
			_processedItems, false);
		mapPosWidth.insert(std::pair<int, std::set<int>>(_processedItems[idx]->idx, possiblePositionsWidth));
		mapPosHeight.insert(std::pair<int, std::set<int>>(_processedItems[idx]->idx, possiblePositionsHeight));
	}
	auto lb = solve(_processedItems, mapPosWidth, mapPosHeight,false) + _processedH;
	if (std::floor(lb) - lb < BLEU::tolerance) lb = std::floor(lb);
	else lb = std::ceil(lb);
	return lb;
	
}




/*
if it should be rotated then rotate the instance, if not do nothing
*/

const bool StripPacking::BLEU::ifRotateInstance(const std::vector<item*>& t_items,const int t_binHeight, const int t_binWidth) const
{
	int sumH = 0;
	int sumW = 0;
	for (size_t idx = 0; idx < t_items.size(); ++idx)
	{
		auto possiblePositionsWidth = computeFX(t_binWidth - t_items[idx]->width, idx,
			t_items, true);
		auto possiblePositionsHeight = computeFX(t_binHeight - t_items[idx]->height, idx,
			t_items, false);
		sumH  += possiblePositionsHeight.size();
		sumW += possiblePositionsWidth.size();
	}
	return (sumH < sumW);
}


void StripPacking::BLEU::rotateInstance(std::vector<item*>& t_Items, int& t_binWidth, int& t_binHeight) const
{
	for (const auto& it : t_Items)
	{
		int tmpW = it->width;
		it->width = it->height;
		it->height = tmpW;
	}
	int tmpWidth = t_binWidth;
	t_binWidth = t_binHeight;
	t_binHeight = tmpWidth;
}



const double StripPacking::BLEU::DualFeasibleFunction1(const int t_alpha, const int t_width) const
{
	double tmp = (t_alpha + 1.0)*((double(t_width) / _processedW));
	double intPart;
	if (std::modf(tmp, &intPart) < BLEU::tolerance) // if tmp is an integer
		return t_width;
	else
		return (std::floor(tmp)*(_processedW / double(t_alpha)));
}


const double StripPacking::BLEU::DualFeasibleFunction2(const int t_alpha, const int t_width) const
{
	if (t_width > _processedW - t_alpha)
		return _processedW;
	else
	{
		if (t_alpha <= t_width && t_width <= _processedW - t_alpha)
			return t_width;
		else return 0;
	}
}

const	 double StripPacking::BLEU::DualFeasibleFunction3(const int t_alpha, const int t_width) const
{
	if (t_width > (_processedW / 2.0))
		return 2 * (std::floor(_processedW / double(t_alpha)) - std::floor((_processedW - t_width) / double(t_alpha)));
	if (std::abs(t_width - (_processedW / 2.0)) < StripPacking::BLEU::tolerance)
		return std::floor(_processedW / double(t_alpha));
	if (t_width < (_processedW / 2.0))
		return 2 * std::floor(t_width / double(t_alpha));
}


void StripPacking::BLEU::dumpSolution(const char* file_name) const
{
	std::ofstream ofs(file_name);
	std::ostringstream ss;
	for (size_t i = 0; i < _processedItems.size(); ++i)
	{
		auto tmp = _processedItems[i];
		ss.str("");
		ss.clear();
		ss << tmp->idx << "," << _finalSolution[tmp->idxHelper].x << "," << _finalSolution[tmp->idxHelper].y << "\n";
		ofs << ss.str();
	}
	ofs.close();
	ofs.open("Rectangle.output");
	for (size_t i = 0; i < _allItems.size(); ++i)
	{
		ss.str("");
		ss.clear();
		ss << _allItems[i]->idx << "," << _allItems[i]->width << "," << _allItems[i]->height << "\n";
		ofs << ss.str();
	}
	ofs.close();
}

void StripPacking::BLEU::dumpSolution(const char* file_name,const std::vector<const item*>& t_Items,
	const std::vector<coordinate>& t_Solution) const
{
	std::ofstream ofs(file_name);
	std::ostringstream ss;
	for (size_t i = 0; i < t_Items.size(); ++i)
	{
		auto tmp = t_Items[i];
		ss.str("");
		ss.clear();
		ss << tmp->idx << "," << t_Solution[tmp->idxHelper].x << "," << t_Solution[tmp->idxHelper].y << "\n";
		ofs << ss.str();
	}
	ofs.close();
	ofs.open("Rectangle.output");
	for (size_t i = 0; i < t_Items.size(); ++i)
	{
		ss.str("");
		ss.clear();
		ss << t_Items[i]->idx << "," << t_Items[i]->width << "," << t_Items[i]->height << "\n";
		ofs << ss.str();
	}
	ofs.close();
}

void StripPacking::BLEU::dumpSolution(const char* file_name, const std::vector<item*>& t_Items,
	const std::vector<coordinate>& t_Solution) const
{
	std::ofstream ofs(file_name);
	std::ostringstream ss;
	for (size_t i = 0; i < t_Items.size(); ++i)
	{
		auto tmp = t_Items[i];
		ss.str("");
		ss.clear();
		ss << tmp->idx << "," << t_Solution[tmp->idxHelper].x << "," << t_Solution[tmp->idxHelper].y << "\n";
		ofs << ss.str();
	}
	ofs.close();
	ofs.open("Rectangle.output");
	for (size_t i = 0; i < t_Items.size(); ++i)
	{
		ss.str("");
		ss.clear();
		ss << t_Items[i]->idx << "," << t_Items[i]->width << "," << t_Items[i]->height << "\n";
		ofs << ss.str();
	}
	ofs.close();
}
// helper functions
const std::set<int> StripPacking::BLEU::getItemsByCol(const int t_Col, 
	const std::vector<item*>& t_Items, const std::vector<coordinate>& t_Cords) const
{
	std::set<int> result;
	for (const auto& it : t_Items)
	{
		if (t_Cords[it->idxHelper].x == -1) continue;
		if (t_Cords[it->idxHelper].x <= t_Col && t_Col <= t_Cords[it->idxHelper].x + it->width - 1)
			result.insert(it->idxHelper);
	}
	return result;
}