#ifndef NODE_H_
#define NODE_H_

#include "../TypeDef.h"
#include "grid_def.h"

#include <math.h>

namespace Larus {
namespace Grid {

enum NodeIdx {
	//=========================
	//   y
	//   |
	//   ---------------
	//   |  PM  |  PP  |
	//   |  2   |  3   |
	//   |  WN  |  NE  |
	//   ---------------
	//   |  SW  |  SE  |
	//   |  0   |  1   |
	//   |  MM  |  MP  |
	//   ------------------->x
	//
	//   ---------------  ---------------
	//   |  MPM |  MPP |  |  PPM |  PPP |
	//   |  2   |  3   |  |  6   |  7   |
	//   |  WNB |  NEB |  |  WNF |  NEF |
	//   ---------------  ---------------
	//   |  SWB |  SEB |  |  SWP |  SEP |
	//   |  0   |  1   |  |  4   |  5   |
	//   |  MMM |  MMP |  |  PMM |  PMP |
	//   ---------------  ---------------
	//=========================

	//2D
	_MM_ = 0,
	_MP_ = 1,
	_PM_ = 2,
	_PP_ = 3,
	//3D
	_MMM_ = 0,
	_MMP_ = 1,
	_MPM_ = 2,
	_MPP_ = 3,
	_PMM_ = 4,
	_PMP_ = 5,
	_PPM_ = 6,
	_PPP_ = 7,
};

inline bool is_x_p(size_t i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 6) == 7;
}
inline bool is_x_m(size_t i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 6) == 6;
}
inline bool is_y_p(size_t i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 5) == 7;
}
inline bool is_y_m(size_t i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 5) == 5;
}
inline bool is_z_p(size_t i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 3) == 7;
}
inline bool is_z_m(size_t i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 3) == 3;
}

template<typename COO_VALUE, typename VALUE, int DIM>
class Node {
public:
	static const size_t Dim = DIM;
	static const size_t NumFaces = DIM + DIM;
	static const size_t NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const size_t NumNeighbors = NumFaces;
	static const size_t NumChildren = NumVertexes;

	typedef COO_VALUE coo_value_t;
	typedef VALUE value_t;
	typedef Node<COO_VALUE, VALUE, DIM> self;
	typedef Node<COO_VALUE, VALUE, DIM>* pself;
	typedef Cell<COO_VALUE, Dim> cell_t;
	typedef cell_t* pcell;
	typedef Data<VALUE, Dim> data_t;
	typedef data_t* pdata;
	typedef self node;
	typedef self* pnode;
	typedef void (*pfunction)(pnode, utPointer);
	typedef void (*pfunction_conditional)(arrayList&, pnode, utPointer);

protected:
	int _node_type;
	size_t _level;
	size_t _root_idx;
	size_t _path;
public:
	pnode father;
	pnode child[NumChildren];
	pnode neighbor[NumNeighbors];
	pcell cell;
	pdata data;

protected:
	int _height(const pnode Current) const {
		if (Current == NULL_PTR) {
			return 0;
		}
		if (!Current->has_child()) {
			return 0;
		} else {
			arrayListV<size_t> arrh(NumChildren);
			for (size_t i = 0; i < this->NumChildren; ++i) {
				arrh = _height(Current->child[i]);
			}
			return 1 + arrh.max();
		}
	}
	void _traversal_conditional(pnode pn, pfunction_conditional pfun_con,
			pfunction pfun, utPointer utp) {
		if (pn == NULL_PTR) {
			return;
		} else {
			(*pfun)(pn, utp);
			if (pn->has_child()) {
				arrayList avt(NumChildren);
				pfun_con(avt, pn, utp);
				for (int i = 0; i < NumChildren; i++) {
					pnode c = pn->child[i];
					if (c != NULL_PTR && avt[i] == 1) {
						_traversal_conditional(c, pfun_con, pfun, utp);
					}
				}
			}
		}
	}
	void _traversal(pnode pn, pfunction pfun, utPointer utp) {
		if (pn == NULL_PTR) {
			return;
		} else {
			(*pfun)(pn, utp);
			if (pn->has_child()) {
				for (int i = 0; i < NumChildren; i++) {
					pnode c = pn->child[i];
					if (c != NULL_PTR) {
						_traversal(c, pfun, utp);
					}
				}
			}
		}
	}
public:
	/*
	 *  constructor
	 */
	Node(pnode f, int nt, size_t level, size_t root_idx, size_t path,	//
			const value_t& x, const value_t& dhx, //
			const value_t& y = 0.0, const value_t& dhy = 0.0, //
			const value_t& z = 0.0, const value_t& dhz = 0.0) {
		_node_type = nt;
		_level = level;
		cell = new cell_t(x, dhx, y, dhy, z, dhz);
		father = f;
		_root_idx = root_idx;
		_path = path;

		data = NULL_PTR;
		for (int i = 0; i < this->NumChildren; i++) {
			child[i] = NULL_PTR;
		}
		for (int i = 0; i < this->NumNeighbors; i++) {
			neighbor[i] = NULL_PTR;
		}
	}
	/*
	 *  delete
	 */
protected:
	/*
	 *  before using this function, making sure that this node is a leaf
	 */
	void _delete_leaf() {
		pnode f = this->father;
		if (f != NULL_PTR) {
			f->child[get_idx()] = NULL_PTR;
		}
		delete cell;
		if (data != NULL_PTR) {
			delete data;
		}
	}
	void _delete(pnode pn) {
		if (pn == NULL_PTR) {
			return;
		}
		if (pn->has_child()) {
			for (int i = 0; i < NumChildren; i++) {
				pnode ch = pn->child[i];
				if (ch != NULL_PTR) {
					_delete(ch);
				}
			}
		} else { // is leaf
			pn->_delete_leaf();
		}
	}
public:
	~Node() {
		_delete(this);
	}
	/*
	 * type
	 */
	inline int get_type() const {
		return _node_type;
	}
	inline void set_type(int type) {
		_node_type = type;
	}

	inline size_t get_level() const {
		return _level;
	}
	inline size_t get_idx() const {
		return (_path >> int(pow(Dim, _level))) & (NumVertexes - 1);
	}
	inline size_t get_path() const {
		return _path;
	}
	inline size_t get_root_idx() const {
		return _root_idx;
	}
	inline size_t height() const {
		return this->_height(this);
	}
	/*
	 *  child
	 */
	inline bool has_child() const {
		for (size_t i = 0; i < this->NumChildren; ++i) {
			if (this->child[i] != NULL_PTR) {
				return true;
			}
		}
		return false;
	}
	inline bool is_leaf() const {
		return !has_child();
	}
	inline bool has_child(size_t idx) const {
		return this->child[idx] != NULL_PTR;
	}
	inline bool is_root() const {
		if (this->father == NULL_PTR) {
			return true;
		} else {
			return false;
		}
	}

	inline bool is_full_child() const {
		bool res = this->child[0] != NULL_PTR;
		for (size_t i = 1; i < this->NumChildren; ++i) {
			res = res && (this->child[i] != NULL_PTR);
		}
		return res;
	}
	inline size_t count_child() const {
		size_t res = 0;
		for (size_t i = 0; i < this->NumChildren; ++i) {
			res += (this->child[i] != NULL_PTR) ? 1 : 0;
		}
		return res;
	}
	/*
	 *  new
	 */
	void new_full_child() {
		if (!has_child()) {
			size_t ltmp = _level + 1;
			value_t nhdx = this->cell->get_hd(_X_) * 0.5;
			value_t nhdy = this->cell->get_hd(_Y_) * 0.5;
			value_t nhdz = this->cell->get_hd(_Z_) * 0.5;
			value_t cx = this->cell->get(_C_, _X_);
			value_t cy = this->cell->get(_C_, _Y_);
			value_t cz = this->cell->get(_C_, _Z_);
			for (size_t i = 0; i < this->NumChildren; ++i) {
				pnode f = this;
				int nt = 1;
				size_t l = ltmp;
				size_t ridx = _root_idx;
				size_t npath = (i << int((pow(Dim, l)))) + _path;
				this->child[i] = new node( //
						f, nt, l, ridx, npath, //
						cx + (is_x_p(i) ? nhdx : -nhdx), nhdx, //
						cy + (is_y_p(i) ? nhdx : -nhdx), nhdy, //
						cz + (is_z_p(i) ? nhdx : -nhdx), nhdz);
			}
		}
	}
	/*
	 *  neighbor find
	 */
	inline bool is_adjacent(const Direction& d) const {
		// Direction on x y or z
		size_t hi = d >> 3;
		return (hi & get_idx()) ^ (hi & d) == 0;
	}
	inline size_t reflect(const Direction& d) const {
		// Direction on x y or z
		return get_idx() ^ (d >> 3);
	}
	inline bool has_diagonal_sibling(const Direction& d) const {
		return (get_idx() ^ (d >> 3)) == (d & 7);
	}
	inline bool is_out_corner(const Direction& d) const {
		return get_idx() == (d & 7);
	}
	inline size_t out_common_direction(const Direction& d) const {
		// return direction on x y or z
		size_t hi = d >> 3;
		size_t low = d & 7;
		return (((low ^ get_idx()) ^ hi) << 3) + low;
	}

};

	/*
	 *  functions out of class
	 */
template <template<typename COO_VALUE, typename VALUE, int DIM>
int get_data_idx(const Node<COO_VALUE, VALUE, DIM>* pn){
	ASSERT(pn!=NULL_PTR);
	return pn->data->get_idx();
}
}//
}
#endif /* NODE_H_ */
