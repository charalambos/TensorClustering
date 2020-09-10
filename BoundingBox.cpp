////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
//////////////////////////////////////////////////////////////////////////////////// 

#ifndef __BOUNDING_BOX_CPP__
#define __BOUNDING_BOX_CPP__

#include "BoundingBox.h"

BoundingBox::BoundingBox()	{
	min_pt = Vector3f(FLT_MAX, FLT_MAX, FLT_MAX);
	max_pt = Vector3f(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	onb = 0x00;
}

BoundingBox::BoundingBox(Vector3f const &_min_pt, Vector3f const &_max_pt, ONB *_onb)	{
	min_pt = _min_pt;
	max_pt = _max_pt;
	onb = _onb;
}

BoundingBox::~BoundingBox()	{
	if (onb)	delete onb;
}


Vector3f BoundingBox::getMinPt() const	{
	return min_pt;
}

Vector3f BoundingBox::getMaxPt() const	{
	return max_pt;
}

ONB *BoundingBox::getONB() const	{
	return onb;
}

void BoundingBox::setMinPt(Vector3f const _min_pt)	{
	min_pt = _min_pt;
	return;
}

void BoundingBox::setMaxPt(Vector3f const _max_pt)	{
	max_pt = _max_pt;
	return;
}

void BoundingBox::setONB(ONB *_onb)	{
	if (onb)	delete onb;
	onb = _onb;
	return;
}


#endif
