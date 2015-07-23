/*
 * ts_segment.h
 *
 *  Created on: May 31, 2015
 *      Author: zhou
 */

#ifndef TS_SEGMENT_H_
#define TS_SEGMENT_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"

namespace LarusTS {

template<class TYPE, std::size_t DIM> class Segment;
template<class TYPE, std::size_t DIM> class Surface;
template<class TYPE, std::size_t DIM> class Edge;

template<class TYPE, std::size_t DIM>
class Segment {
public:
	typedef Segment<TYPE, DIM> self_class;
	typedef std::size_t size_type;
	typedef Point<TYPE, DIM> Poi;
	typedef Point<TYPE, DIM>* pPoi;
	typedef Segment<TYPE, DIM> Seg;
	typedef Segment<TYPE, DIM>* pSeg;
	typedef Vertex<TYPE, DIM> Ver;
	typedef Vertex<TYPE, DIM>* pVer;
	typedef std::list<pSeg> list_pSeg;
	typedef std::list<pVer> list_pVer;
public:
	pVer v1;
	pVer v2;
public:
	Segment() {
		v1 = nullptr;
		v2 = nullptr;
	}
	Segment(pVer a, pVer b) {
		assert(a != nullptr);
		assert(b != nullptr);
		v1 = a;
		v2 = b;
		v1->segments.push_back(this);
		v2->segments.push_back(this);
	}

};

//function out of class
/**
 * gts_segments_are_intersecting:
 * @s1: a #GtsSegment.
 * @s2: a #GtsSegment.
 *
 * Returns: %GTS_IN if @s1 and @s2 are intersecting, %GTS_ON if one of the
 * endpoints of @s1 (resp. @s2) lies on @s2 (resp. @s1), %GTS_OUT otherwise.
 */
template<class TYPE, std::size_t DIM>
Intersect segments_are_intersecting( //
		const Segment<TYPE, DIM>& s1, //
		const Segment<TYPE, DIM>& s2) {
	typename Segment<TYPE, DIM>::Poi * p1, *p2, *p3, *p4;
	Float d1, d2, d3, d4;

	p1 = s1.v1;
	p2 = s1.v2;
	p3 = s2.v1;
	p4 = s2.v2;
	d1 = point_orientation(p1, p2, p3);
	d2 = point_orientation(p1, p2, p4);
	if ((d1 > 0.0 && d2 > 0.0) || (d1 < 0.0 && d2 < 0.0))
		return TS_OUT;
	d3 = point_orientation(p3, p4, p1);
	d4 = point_orientation(p3, p4, p2);
	if ((d3 > 0.0 && d4 > 0.0) || (d3 < 0.0 && d4 < 0.0))
		return TS_OUT;
	if (d1 == 0.0 || d2 == 0.0 || d3 == 0.0 || d4 == 0.0)
		return TS_ON;
	return TS_IN;
}

}

#endif /* TS_TS_SEGMENT_H_ */
