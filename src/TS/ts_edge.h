/*
 * ts_edge.h
 *
 *  Created on: May 31, 2015
 *      Author: zhou
 */

#ifndef TS_EDGE_H_
#define TS_EDGE_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"

namespace LarusTS {

template<class TYPE, std::size_t DIM> class Triangle;
template<class TYPE, std::size_t DIM> class Surface;

template<class TYPE, std::size_t DIM>
class Edge: public Segment<TYPE, DIM> {
public:
	typedef Segment<TYPE, DIM> base_class;
	typedef Edge<TYPE, DIM> self_class;
	typedef std::size_t size_type;
	typedef Point<TYPE, DIM> Poi;
	typedef Point<TYPE, DIM>* pPoi;
	typedef Segment<TYPE, DIM> Seg;
	typedef Segment<TYPE, DIM>* pSeg;
	typedef Vertex<TYPE, DIM> Ver;
	typedef Vertex<TYPE, DIM>* pVer;
	typedef Triangle<TYPE, DIM> Tri;
	typedef Triangle<TYPE, DIM>* pTri;
	typedef std::list<pSeg> list_pSeg;
	typedef std::list<pVer> list_pVer;
	typedef std::list<pTri> list_pTri;
public:
	list_pTri triangles;

	Edge(pVer a, pVer b) :
			base_class(a, b) {
	}

};

}

#endif /* TS_TS_EDGE_H_ */
