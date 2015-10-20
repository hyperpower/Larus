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
	typedef Node<COO_VALUE, VALUE, DIM> Self;
	typedef Node<COO_VALUE, VALUE, DIM>* pSelf;
	typedef Cell<COO_VALUE, Dim> Cell_;
	typedef Cell_* pCell;
	typedef Data<VALUE, Dim> Data_;
	typedef Data_* pData;
	typedef Self Node_;
	typedef Self* pNode;
	typedef void (*pFun)(pNode, utPointer);
	typedef void (*pFun_Conditional)(arrayList&, pNode, utPointer);

protected:
	typedef size_t st;
	typedef COO_VALUE cvt;
	typedef VALUE vt;
	//
	int _node_type;
	st _level;
	st _root_idx;
	st _path;
public:
	pNode father;
	pNode child[NumChildren];
	pNode neighbor[NumNeighbors];
	pCell cell;
	pData data;

protected:
	int _height(const pNode Current) const {
		if (Current == NULL_PTR) {
			return 0;
		}
		if (!Current->HasChild()) {
			return 0;
		} else {
			arrayListV<st> arrh(NumChildren);
			for (st i = 0; i < this->NumChildren; ++i) {
				arrh = _height(Current->child[i]);
			}
			return 1 + arrh.max();
		}
	}
	void _traversal_conditional(pNode pn, pFun_Conditional pfun_con,
			pFun pfun, utPointer utp) {
		if (pn == NULL_PTR) {
			return;
		} else {
			(*pfun)(pn, utp);
			if (pn->HasChild()) {
				arrayList avt(NumChildren);
				pfun_con(avt, pn, utp);
				for (int i = 0; i < NumChildren; i++) {
					pNode c = pn->child[i];
					if (c != NULL_PTR && avt[i] == 1) {
						_traversal_conditional(c, pfun_con, pfun, utp);
					}
				}
			}
		}
	}
	void _traversal(pNode pn, pFun pfun, utPointer utp) {
		if (pn == NULL_PTR) {
			return;
		} else {
			(*pfun)(pn, utp);
			if (pn->HasChild()) {
				for (int i = 0; i < NumChildren; i++) {
					pNode c = pn->child[i];
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
	Node(pNode f, int nt, st level, st root_idx, st path,	//
			const vt& x, const vt& dhx, //
			const vt& y = 0.0, const vt& dhy = 0.0, //
			const vt& z = 0.0, const vt& dhz = 0.0) {
		_node_type = nt;
		_level = level;
		cell = new Cell_(x, dhx, y, dhy, z, dhz);
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
	void _DeleteLeaf() {
		pNode f = this->father;
		if (f != NULL_PTR) {
			f->child[GetIdx()] = NULL_PTR;
		}
		delete cell;
		if (data != NULL_PTR) {
			delete data;
		}
	}
	void _Delete(pNode pn) {
		if (pn == NULL_PTR) {
			return;
		}
		if (pn->HasChild()) {
			for (int i = 0; i < NumChildren; i++) {
				pNode ch = pn->child[i];
				if (ch != NULL_PTR) {
					_Delete(ch);
				}
			}
		} else { // is leaf
			pn->_DeleteLeaf();
		}
	}
public:
	~Node() {
		_Delete(this);
	}
	/*
	 * type
	 */
	inline int GetType() const {
		return _node_type;
	}
	inline void SetType(int type) {
		_node_type = type;
	}

	inline st GetLevel() const {
		return _level;
	}
	inline st GetIdx() const {
		return (_path >> int(pow(Dim, _level))) & (NumVertexes - 1);
	}
	inline st GetPath() const {
		return _path;
	}
	inline st GetRootIdx() const {
		return _root_idx;
	}
	inline st Height() const {
		return this->_height(this);
	}
	/*
	 *  child
	 */
	inline bool HasChild() const {
		for (st i = 0; i < this->NumChildren; ++i) {
			if (this->child[i] != NULL_PTR) {
				return true;
			}
		}
		return false;
	}
	inline bool IsLeaf() const {
		return !HasChild();
	}
	inline bool HasChild(st idx) const {
		return this->child[idx] != NULL_PTR;
	}
	inline bool IsRoot() const {
		if (this->father == NULL_PTR) {
			return true;
		} else {
			return false;
		}
	}

	inline bool IsFullChild() const {
		bool res = this->child[0] != NULL_PTR;
		for (st i = 1; i < this->NumChildren; ++i) {
			res = res && (this->child[i] != NULL_PTR);
		}
		return res;
	}
	inline st CountChild() const {
		st res = 0;
		for (st i = 0; i < this->NumChildren; ++i) {
			res += (this->child[i] != NULL_PTR) ? 1 : 0;
		}
		return res;
	}
	/*
	 *  new
	 */
	void NewFullChild() {
		if (!HasChild()) {
			st ltmp = _level + 1;
			vt nhdx = this->cell->get_hd(_X_) * 0.5;
			vt nhdy = this->cell->get_hd(_Y_) * 0.5;
			vt nhdz = this->cell->get_hd(_Z_) * 0.5;
			vt cx = this->cell->get(_C_, _X_);
			vt cy = this->cell->get(_C_, _Y_);
			vt cz = this->cell->get(_C_, _Z_);
			for (st i = 0; i < this->NumChildren; ++i) {
				pNode f = this;
				int nt = 1;
				st l = ltmp;
				st ridx = _root_idx;
				st npath = (i << int((pow(Dim, l)))) + _path;
				this->child[i] = new Node_( //
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
	inline bool IsAdjacent(const Direction& d) const {
		// Direction on x y or z
		st hi = d >> 3;
		return (hi & GetIdx()) ^ (hi & d) == 0;
	}
	inline st Reflect(const Direction& d) const {
		// Direction on x y or z
		return GetIdx() ^ (d >> 3);
	}
	inline bool HasDiagonalSibling(const Direction& d) const {
		return (GetIdx() ^ (d >> 3)) == (d & 7);
	}
	inline bool IsOutCorner(const Direction& d) const {
		return GetIdx() == (d & 7);
	}
	inline st OutCommonDirection(const Direction& d) const {
		// return direction on x y or z
		st hi = d >> 3;
		st low = d & 7;
		return (((low ^ GetIdx()) ^ hi) << 3) + low;
	}

};

	/*
	 *  functions out of class
	 */
template <typename COO_VALUE, typename VALUE, int DIM>
int GetDataIdx(const Node<COO_VALUE, VALUE, DIM>* pn){
	ASSERT(pn!=NULL_PTR);
	return pn->data->get_idx();
}
}//
}
#endif /* NODE_H_ */
