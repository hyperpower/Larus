#ifndef GRID_H_
#define GRID_H_

#include "../TypeDef.h"
#include "grid_def.h"
#include "node.h"
#include "cell.h"

#include "../Algebra/Space.h"

namespace Larus {
namespace Grid {
template<typename COO_VALUE, typename VALUE, int DIM>
class Grid {
public:
	static const size_t Dim = DIM;
	static const size_t NumFaces = DIM + DIM;
	static const size_t NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const size_t NumNeighbors = NumFaces;

	typedef COO_VALUE coo_value_t;
	typedef VALUE value_t;
	typedef Grid<COO_VALUE, VALUE, DIM> self;
	typedef Grid<COO_VALUE, VALUE, DIM>* pself;
	typedef Cell<COO_VALUE, Dim> cell_t;
	typedef cell_t* pcell;
	typedef Data<VALUE, Dim> data_t;
	typedef data_t* pdata;
	typedef Node<COO_VALUE, VALUE, DIM> node;
	typedef Node<COO_VALUE, VALUE, DIM>* pnode;
	typedef void (*pfunction)(pnode, utPointer);
	typedef void (*pfunction_conditional)(arrayList&, pnode, utPointer);
	/*
	 *  data
	 */
	SpaceT<pnode, Dim> grid;
	/*
	 *  constructor
	 */
	Grid() {

	}
	Grid(size_t ni, coo_value_t ox, coo_value_t dx, //
			size_t nj = 0, coo_value_t oy = 0, coo_value_t dy = 0, //
			size_t nk = 0, coo_value_t oz = 0, coo_value_t dz = 0) {
		grid.recontruct(ni, nj, nk);
		for (int i = 0; i < ni; i++) {
			for (int j = 0; j < (Dim >= 2) ? nj : 1; j++) {
				for (int k = 0; k < (Dim == 3) ? nk : 1; k++) {
					// new
					grid.at_1d(i) =  //
							new node(
							NULL_PTR, // father
									0, // type
									0, //level
									0, //node idx
									ox + (i + 0.5) * dx, 0.5 * dx, //x
									oy + (i + 0.5) * dy, 0.5 * dy, //y
									oz + (i + 0.5) * dz, 0.5 * dz); //z
				}
			}
		}
	}
protected:
	void _delete() {
		for (int i = 0; i < grid.size(); i++) {
			if (grid.at_1d(i) != NULL_PTR) {
				delete grid.at_1d(i);
			}
		}
	}
public:
	~Grid() {
		_delete();
	}
	/*
	 *  size
	 */
	inline size_t i_len() const {
		return grid.iLen();
	}
	inline size_t j_len() const {
		return grid.jLen();
	}
	inline size_t k_len() const {
		return (Dim < 3) ? 0 : grid.kLen();
	}
	inline bool is_empty() const {
		if (grid.size() <= 0) {
			return true;
		} else {
			return false;
		}
	}
	size_t get_dim() const {
		return Dim;
	}
	size_t size() const {
		return grid.size();
	}

};
}
}
#endif
