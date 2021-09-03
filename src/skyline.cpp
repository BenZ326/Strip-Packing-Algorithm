#include "skyline.h"
#include "heuristic.h"
// select a skyline for an item to place according to the given mode
StripPacking::Skyline* StripPacking::selectSkyline(const StripPacking::Skyline* t_head, 
	const StripPacking::skylineSelectionMode & t_mode)
{
	switch (t_mode)
	{
	case StripPacking::skylineSelectionMode::leftBottom:
	{
		Skyline* cur = t_head->next;
		Skyline* res = nullptr;
		int min_height = 999999;
		while (cur != nullptr)
		{
			if (cur->corY < min_height) {
				min_height = cur->corY;
				res = cur;
			}
			cur = cur->next;
		}
		return res;
	}
	case StripPacking::skylineSelectionMode::bestFit:
	{
		Skyline* cur = t_head->next;
		Skyline* res = nullptr;
		int min_height = 999999;
		while (cur != nullptr)
		{
			if (cur->corY < min_height) {
				min_height = cur->corY;
				res = cur;
			}
			cur = cur->next;
		}
		return res;
	}
	default:
		break;
	}
}

// add an item over the selected skyline, see the relaxationMode
void StripPacking::addItemOverSkyline(StripPacking::Skyline* t_skyline,
	const StripPacking::item* t_item)
{
	// update the placement of the solution
	Heuristic::solutions[t_item->idx].x = t_skyline->corX;
	Heuristic::solutions[t_item->idx].y = t_skyline->corY;
	// update the skyline
	const int scenario = t_item->width < t_skyline->length ? 1 : (t_item->width == t_skyline->length ? 0 : 2);
	switch (scenario)
	{
	case 0:
	{
		t_skyline->corY += t_item->height;
		detectAndMergeSkylines(t_skyline);
		break;

	}
	case 1:
	{
		Skyline* new_skyline = new Skyline();
		new_skyline->corX = t_skyline->corX + t_item->width;
		new_skyline->corY = t_skyline->corY;
		t_skyline->corY += t_item->height;
		new_skyline->length = t_skyline->length - t_item->width;
		t_skyline->length = t_item->width;
		new_skyline->next = t_skyline->next;
		new_skyline->next->prev = new_skyline;
		t_skyline->next = new_skyline;
		new_skyline->prev = t_skyline;
		detectAndMergeSkylines(t_skyline);
		break;
	}
	case 2:
	{
		break;
	}
	default:
		break;
	}

}


void StripPacking::detectAndMergeSkylines(StripPacking::Skyline* t_skyline)
{
	// leftwards
	auto cur = t_skyline->prev;
	while (cur != nullptr && cur->corY == t_skyline->corY)
	{
		cur->length += t_skyline->length;
		removeSkyline(t_skyline);
		t_skyline = cur;
		cur = cur->prev;
	}
	// rightwards
	cur = t_skyline->next;
	while (cur != nullptr && cur->corY == t_skyline->corY)
	{
		t_skyline->length += cur->length;
		removeSkyline(cur);
		cur = t_skyline->next;
	}
}

void StripPacking::removeSkyline(Skyline* t_skyline)
{
	auto pre = t_skyline->prev;
	auto nxt = t_skyline->next;
	delete t_skyline;
	pre->next = nxt;
	nxt->prev = pre;
}


void StripPacking::liftSkyline(StripPacking::Skyline* t_skyline)				// lift a skyline to the adjacent  skyline that has lower height
{
	// left
	if (t_skyline->prev->corY < t_skyline->next->corY)
	{
		t_skyline->prev->length += t_skyline->length;
		removeSkyline(t_skyline);
		return;
	}
	if (t_skyline->prev->corY == t_skyline->next->corY )
	{
		if (t_skyline->prev->corY == 999999)
		{
			t_skyline->corY = 999999;
		}
		else
		{
			t_skyline->prev->length += t_skyline->length;
			auto pre = t_skyline->prev;
			removeSkyline(t_skyline);
			detectAndMergeSkylines(pre);
		}
		return;
	}
	if (t_skyline->prev->corY > t_skyline->next->corY)
	{
		t_skyline->next->length += t_skyline->length;
		t_skyline->next->corX = t_skyline->corX;
		removeSkyline(t_skyline);
		return;
	}
}