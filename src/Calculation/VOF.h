/************************
 //  \file   VOF.h
 //  \brief
 // 
 //  \author czhou
 //  \date   27 janv. 2015 
 ***********************/
#ifndef VOF_H_
#define VOF_H_

#include "../TypeDef.h"
#include "CalDef.h"
#include "../Geometry/Plane.h"
#include "../Utility/Pair.h"
#include "../Utility/List.h"
#include "Scalar.h"

namespace Larus {
//Function for 3D================================
//===============================================
int whichcase_3D_8(Float, Float, Float);

Float newAlpha_Forward(const Plane&, Float, Float, Float);

Float calFractionVolume(const Plane&, Float, Float, Float);

void coe_premutation( //
		const Pair<Float, Float>& o1, //
		const Pair<Float, Float>& o2, //
		const Pair<Float, Float>& o3, //
		Pair<Float, Float>& p1, // out
		Pair<Float, Float>& p2, // out
		Pair<Float, Float>& p3); //out

ListT<Point3D> calInterctPoints(const Plane &l, Float c1, Float c2, Float c3);
void transfer_ListP(ListT<Point3D> &list, Float x1, Float x2, Float x3);
void scale_ListP(ListT<Point3D> &list, Float x1, Float x2, Float x3);
//end of 3D function ============================

//Function for 2D================================
//===============================================

Float calAlpha(Float a, Float b, Float C);

Line calLine(Float a, Float b, Float C);

Float calFractionArea(const Line &l, Float c1, Float c2);

Float calFractionArea(const Line &l, Float xo, Float yo,Float c1, Float c2);

Segment2D calInterctPoints(const Line &l, Float c1, Float c2);

void getListSegment( //
		ListT<Segment2D>& lseg, //list Segment
		Forest2D& forest,        //forest
		LarusDef::size_type idc, //data color index
		LarusDef::size_type idx, //data index
		LarusDef::size_type idy, //data index
		LarusDef::size_type ida  //data index
		);

Float calInterfaceLength(const ListT<Segment2D> &sg);

//draw===========================================
void draw_vtu_vof( // draw vof 3D==============
		std::string filename,  //filename
		pOCTree oct,           //tree
		LarusDef::size_type idc, //data color index
		LarusDef::size_type idx, //data index
		LarusDef::size_type idy, //data index
		LarusDef::size_type idz, //data index
		LarusDef::size_type ida  //data index
		);
void draw_vtu_vof( // draw vof 3D==============
		std::string filename,  //filename
		Forest3D& forest,      //tree
		LarusDef::size_type idc, //data color index
		LarusDef::size_type idx, //data index
		LarusDef::size_type idy, //data index
		LarusDef::size_type idz, //data index
		LarusDef::size_type ida  //data index
		);
}

#endif /* VOF_H_ */
