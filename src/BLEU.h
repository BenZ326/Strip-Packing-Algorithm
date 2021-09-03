#pragma once
#include "spp.h"
#include <ilcplex/ilocplex.h>

#include <stack>
class itemPieceWidth;
namespace StripPacking
{

/*
Fully reproduce the paper: Combinatorial Benders' Cuts for the strip packing problem
*/
class BLEU
{
	// global parameters
public:
	static double tolerance;
	static int bigNumber;
	static int BBMaxExplNodesPerPack;				// maximal number of explored nodes for perfect packing (for the BB algorithm)
	static int BBMaxExplNodesNonPerPack;		// maximal number of explored nodes for non-perfect packing (for the BB algorithm)
	static int interestingStatics;
	static int ycheckExplNode;
	static bool nodeLimitFlag;			// if y-check subroutine reaches node limit, it becomes true;
	static algorithmStatus algStatus;
	/*
	Explanation on the nodeLimitFlag and algStatus;
	Before running an algorithm, 
	initialize with:
		nodeLimitFlag = false; algStatus = exact
	in-process:
		as soon as ycheck exceeds node limit, nodeLimitFlag is flipped to True. 
	end:
		if solution status is infeasible and nodeLimitFlag = true, then the algorithm is no longer exact,
		so algStatus = approximate
	*/
public:
	BLEU(const std::vector<const item*>& t_items, const int t_W, const int t_TrialHeight, const int t_timeLimit);
	BLEU(const std::vector<const item*>& t_items, const int t_W, const int t_timeLimit);
	void preprocessing();
	void bounds();
	const solutionStatus evaluate();
	const StripPacking::solutionStatus solvePCC();			// solve the parallel machine scheduling with contiguity constraints
	void dumpSolution(const char* file_name = "Solution.output") const;
	void dumpSolution(const char* file_name, const std::vector<const item*>& t_Items,
		const std::vector<coordinate>& t_Solution) const;
	void dumpSolution(const char* file_name, const std::vector<item*>& t_Items,
		const std::vector<coordinate>& t_Solution) const;
// algorithms
protected:
	//preprocessing and bounds	
protected:
	// 5.1 preprocessing
	void reassignItemsIdx();
	void preprocessingFixItems();
	void preprocessingReduceW();
	void preprocessingModifyItemWidth();
	void preprocessItemHeight(std::vector<item*>& t_items, 
		const int t_binHeight, int& t_binWidth);			// only after the trial height is determined
	// 5.1 preprocess the bounds
	const int LowerBound1() const;			// 5.1 lower bound 1
	const int LowerBound2() const;			// 5.2 lower bound 2
	const int LowerBound3() const;			// 5.3 lower bound 3
	const int LowerBound4() const;			// 5.4 lower bound 4
	const int LowerBound5() const;			// 5.5 lower bound 5
	// cut items into pieces each of which has 1 unit height and the original width
	void cutItemsAlongHeight();
	// dual feasibility function
const	double DualFeasibleFunction1(const int t_alpha, const int t_width) const;
const double DualFeasibleFunction2(const int t_alpha, const int t_width) const ;
const	double DualFeasibleFunction3(const int t_alpha, const int t_width) const;




// rotate instances
void rotateInstance(std::vector<item*>& t_Items, int& t_binWidth, int& t_binHeight) const;
const bool ifRotateInstance(const std::vector<item*>& t_items, const int t_binHeight, const int t_binWidth) const;
/*
The branch and bound algorithms----------------------------------------start
*/

class oneDimensionItem
{
public:
	oneDimensionItem(const int& t_weight, const int& t_value) :weight(t_weight), value(t_value) {}
	const int weight;
	const int value;
};
class BBNode
{
public:
	BBNode(const std::vector<const item*>& t_remainingItems, const int t_Width, const int t_TrialHeight);
	BBNode(const std::vector<const item*>& t_remainingItems, const std::vector<coordinate>& t_Cords,
		const int t_Width, const int t_TrialHeight);
	BBNode(const BBNode& t_BBNode);					// copy constructor
	BBNode(const BBNode& t_BBNode, const item* t_placedItem);			// invoked when making a branch, the t_placedItem is the packed item in this time of branching
	const int trialHeight;
	int leftMostIdx;
	std::vector<int> columnsOccupiedHeight;			// [10,5,3,2] means 10 units of height in the 1st column is occupied and 5 units for the 2nd column...
	std::vector<const item*> remainingItems;			// make heap
	std::vector<const item*> packedItems;
	std::vector<int> maxiItemIdxColumns;				// [3,5,1,7] means among items placed in 1st column, 3 is the largest index of...
	std::vector<coordinate> itemPositions;				// store the final positions of all the items in processedItems (respect the order in processedItems)
};
protected:
	const solutionStatus branchAndBound(const std::vector<const item*>& t_Items, const int t_binWidth,
		const int t_binHeight);
	const solutionStatus branchAndBoundYRelax(const std::vector<const item*>& t_Items, const int t_binWidth,
		const int t_binHeight);
	void makeBranch(const std::unique_ptr<BBNode>& t_currentNode, 
		std::stack<std::unique_ptr<BBNode>>& t_dfstree) const;
	const bool bounding(const std::unique_ptr<BBNode>& t_currentNode) const;
	const bool dynamicCuts(const std::unique_ptr<BBNode>& t_currentNode) const;

/*
The branch and bound algorithms----------------------------------------end
*/


/*
the y-check algorithm---------------------------------------------------start
*/
bool yCheckAlgorithm(const int t_processedW, const int t_TrialHeight, const std::vector<coordinate>& itemPositions,
	const std::vector<const item*> t_processedItems) const;

/*
The enumerate tree described right before section 4
*/
solutionStatus yCheckEnumerationTree(const std::vector<const item*>& t_InterestItems, const std::vector<coordinate>& t_Cords,
	const int t_Height, const int t_Width) const ;
bool yCheckBounding(const std::unique_ptr<BBNode>& t_currentNode) const;
void yCheckMakeBranch(const std::unique_ptr<BBNode>& t_currentNode, 
	std::stack<std::unique_ptr<BBNode>>& t_yEntree) const;
/*
Given column heights, return a Niche which is featured as the start column and the end column
the return value is a two dimensional array, one element of the array is the start column of the niche
while the second is the end column of the niche

*/
std::vector<int> getNiche(const std::vector<int>& t_ColumnHeights) const;

/*
To preprocess all the items considered in the branch and bound tree to reduce the complexity of the further y-check problem 
*/
const std::vector<const item*> preprocess4yCheck(int& t_Width, const std::vector<const item*>& t_InterestItems, std::vector<coordinate>& t_Cords,
	const int t_Height) const;

const std::vector<item*> preprocessedFirst4yCheck(const std::vector<const item*>& t_InterestItems, 
	std::vector<coordinate>& t_Cords, const int t_Width) const;

void  preprocessedSecond4yCheck(std::vector<item*>& t_allItems,
	std::vector<coordinate>& t_Cords, const int t_binWidth) const;

/*
The function would modify the bin width as well as the width for items
*/
void  preprocessedThird4yCheck(std::vector<item*>& t_allItems,
	std::vector<coordinate>& t_Cords, int& t_binWidth) const;
/*
// input: i, items in the left of i, coordinates, start column, right or left (true is left, false is right)
// output: bool, change the items in the left of i
*/
bool mergeItems4yCheck(std::vector<item*>& t_allItems,item* t_i, std::map<int, std::list<item*>>& t_leftItems, 
	std::map<int, std::list<item*>>& t_rightItems
	,std::list<item*>& t_Items, std::vector<coordinate>& t_Cords,
	const bool t_left, const int t_binWidth) const;

const int getFirstColumn(const std::list<item*>& t_lefts, const std::vector<coordinate>& t_Cords) const;
const int getFirstColumn(const int t_prevStartColumn,const std::list<item*>& t_lefts, const std::vector<coordinate>& t_Cords) const;
const int getLastColumn(const std::list<item*>& t_rights, const std::vector<coordinate>& t_Cords) const;
const int getLastColumn(const int t_prevLastColumn, const std::list<item*>& t_rights, 
	const std::vector<coordinate>& t_Cords) const;
const int getMaxWidth(const int t_column,const std::list<item*>& t_Items, const std::vector<coordinate>& t_Cords, const bool t_left) const;
void transferItemsAndCords4YEnumeration(const std::list<item*>& t_OrigItems, const std::vector<coordinate>& t_OrigCords
	,std::vector<const item*>& t_TransferredItems, std::vector<coordinate>& t_TransferredCords) const;

void merging(item* t_i, std::vector<item*>& t_allItems, std::vector<coordinate>& t_Cords,
	const int t_startColumn, const int t_maxWidth, std::list<item*>& t_Items, const bool t_left) const;
const std::vector<std::map<int, std::list<item*>>> getLeftsAndRights(const std::vector<item*>& t_allItems,
	const std::vector<coordinate>& t_Cords, const std::set<int>& t_coveredItems) const;
void mergeItems(item* t_i, std::list<item*>& t_Items, std::map<int, std::list<item*>>& t_allItems, std::vector<coordinate>& t_Cords) const;
void releaseTmpItems(std::vector<const item*>& t_Items) const;
const bool checkSeparable(const std::list<item*>& t_Items, const std::vector<coordinate>& t_Cords, 
	const int t_startColumn, const int t_maxWidth, const bool t_left) const;
/*
the y-check algorithm---------------------------------------------------end
*/

/*
Return:
	A vector of endColumns of each subset (suppose we have W = 15)
	e.g. [5,9] represents three subsets of columns [0,1,2,3,4,5] [6,7,8,9], [10,11,12,13,14]. 
*/
const std::vector<std::vector<int>> verticalCuts(const int t_binHeight, const int t_binWidth, const std::vector<const item*>& t_allItems, 
	const std::vector<coordinate>& t_Cords) const;



//---- helper functions
/* t_Cords respects the idxHelper , return the set of idxHelpers of items  that occupy the column t_Col*/
const std::set<int> getItemsByCol(const int t_Col, const std::vector<item*>& t_Items, const std::vector<coordinate>& t_Cords) const;
private:
	std::vector<const item*> _allItems;
	std::vector<const itemPieceWidth*> _allItemPiecesWidths;
	std::vector<int> _allWidths;
	std::vector<const item*> _processedItems;				// items that are going to be computed
	const int _W;
	int _processedH;
	int _processedW;
	int _bestLowerBound;
	int _trialHeight;			// the current height being tried
	std::vector<coordinate> _finalSolution;
	const bool _evaluatedMode;			// if it is true, then the algorithm starts from a given height, the mission is to determine if the height is feasible
	const int _timeLimit;

};
}

