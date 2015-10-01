#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "../TypeDef.h"
#include "grid_def.h"
#include "node.h"
#include "cell.h"


namespace Larus {
namespace Grid {



struct GhostID {
	int node_idx;  //the global idx of the origin node
	size_t step;   //the steps of ghost node, we can choose multiple ghost node,
				   // usually step = 0
	Direction direction; //The direction only on x, y or z
};

struct GhostID_compare {
	bool operator()(const GhostID& lhs,
			const GhostID& rhs) const {
		if (lhs.node_idx < rhs.node_idx) {
			return true;
		} else if (lhs.node_idx == rhs.node_idx) {
			return int(lhs.direction) < int(rhs.direction);
		} else if (lhs.direction == rhs.direction) {
			return lhs.step < rhs.step;
		} else {
			return false;
		}
	}
};

template<typename COO_VALUE, typename VALUE, int DIM>
class Boundary {
public:
public:
	static const size_t Dim = DIM;
	static const size_t NumFaces = DIM + DIM;
	static const size_t NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const size_t NumNeighbors = NumFaces;

	typedef COO_VALUE coo_value_t;
	typedef VALUE value_t;
	typedef Boundary<COO_VALUE, VALUE, DIM> self;
	typedef const Boundary<COO_VALUE, VALUE, DIM> const_self;
	typedef Boundary<COO_VALUE, VALUE, DIM>* pself;
	typedef const Boundary<COO_VALUE, VALUE, DIM>* const_pself;
	typedef Cell<COO_VALUE, Dim> cell_t;
	typedef cell_t* pcell;
	typedef Data<VALUE, Dim> data_t;
	typedef data_t* pdata;
	typedef Node<COO_VALUE, VALUE, DIM> node;
	typedef Node<COO_VALUE, VALUE, DIM>* pnode;
	typedef void (*pfunction)(pnode, utPointer);
	typedef void (*pfunction_conditional)(arrayList&, pnode, utPointer);

};


}
}

#endif
