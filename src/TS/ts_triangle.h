/*
 * ts_triangle.h
 *
 *  Created on: May 31, 2015
 *      Author: zhou
 */

#ifndef TS_TRIANGLE_H_
#define TS_TRIANGLE_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"
#include "ts_edge.h"

namespace LarusTS {

template<class TYPE, std::size_t DIM> class Surface;

template<class TYPE, std::size_t DIM>
class Triangle {
public:
	typedef Triangle<TYPE, DIM> self_class;
	typedef std::size_t size_type;
	typedef Point<TYPE, DIM> Poi;
	typedef Point<TYPE, DIM>* pPoi;
	typedef Segment<TYPE, DIM> Seg;
	typedef Segment<TYPE, DIM>* pSeg;
	typedef Edge<TYPE, DIM> Edg;
	typedef Edge<TYPE, DIM>* pEdg;
	typedef Vertex<TYPE, DIM> Ver;
	typedef Vertex<TYPE, DIM>* pVer;
	typedef Triangle<TYPE, DIM> Tri;
	typedef Triangle<TYPE, DIM>* pTri;
	typedef std::list<pSeg> list_pSeg;
	typedef std::list<pVer> list_pVer;
	typedef std::list<pTri> list_pTri;
public:
	pEdg e1;
	pEdg e2;
	pEdg e3;
public:
	Triangle(pEdg a, pEdg b, pEdg c){
		e1=a;
		e2=b;
		e3=c;
	}

};

}

#endif /* TS_TS_TRIANGLE_H_ */
