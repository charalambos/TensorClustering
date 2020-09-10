////////////////////////////////////////////////////////////////////////////////////
// Copyright © Charalambos "Charis" Poullis, charalambos@poullis.org    	  //
// This work can only be used under an exclusive license of the author.           //
//////////////////////////////////////////////////////////////////////////////////// 


#ifndef __BOUNDING_BOX_H__
#define __BOUNDING_BOX_H__

/**		Bounding Box
* This class defines a bounding box and related functions
* 
*/


#include <float.h>
#include "ONB.h"


class BoundingBox	{
	public: 
		///Constructor
		BoundingBox();
		///Constructor
		BoundingBox(Vector3f const &_min_pt, Vector3f const &_max_pt, ONB *_onb = 0x00);
		///Destructor
		~BoundingBox();

		///Get min point
		Vector3f getMinPt() const;

		///Get max point
		Vector3f getMaxPt() const;

		///Get ONB
		ONB *getONB()	const;

		///Set min point
		void setMinPt(Vector3f const _min_pt);

		///Set max point
		void setMaxPt(Vector3f const _max_pt);

		///Set the orientation
		void setONB(ONB *_onb);

	protected:
		///The minimum point
		Vector3f min_pt;
		///The maximum point
		Vector3f max_pt;
		///The orthonormal basis of the bbox
		ONB *onb;
};


#endif
