/************************
 //  \file   Forest.h
 //  \brief
 // 
 //  \author czhou
 //  \date   4 févr. 2015 
 ***********************/
#ifndef FOREST_H_
#define FOREST_H_

#include "../TypeDef.h"

#include "SPTree.h"
#include "../IO/IO.h"
#include "../IO/IO_vtk.h"
#include "../IO/Gnuplot.h"
#include "../Algebra/Matrix.h"
#include "../Algebra/Space.h"
#include "Dimension.h"

namespace Larus {

template<class Dimension, class Cell, class Data, LarusDef::size_type Dim,
		class _Ref, class _Ptr>
class _SPNode_iterator_global;

template<class Dimension, class Cell, class Data, LarusDef::size_type Dim,
		class _Ref, class _Ptr>
class _SPNodeFace_iterator_global;

template<class Dimension, class Cell, class Data, LarusDef::size_type Dim,
		class _Ref, class _Ptr>
class _SPNodeVertex_iterator_global;

static const int ATT_DISABLE = 0;
static const int ATT_ENABLE = 1;

template<typename DIM_TYPE>
class Forest {
public:
	typedef typename DIM_TYPE::Tree Tree;
	typedef typename DIM_TYPE::pTree pTree;
	typedef typename DIM_TYPE::Node Node;
	typedef typename DIM_TYPE::pNode pNode;
	typedef typename DIM_TYPE::Node::Cell_type Cell;
	typedef typename DIM_TYPE::Node::Data_type Data;
	typedef typename DIM_TYPE::size_type size_type;
	static const size_type Dim = DIM_TYPE::DIM;

	typedef SPNodeFace<Cell, Data, Dim> Face;
	typedef SPNodeFace<Cell, Data, Dim>* pFace;

	typedef SPNodeVertex<Cell, Data, Dim> Vertex;
	typedef SPNodeVertex<Cell, Data, Dim>* pVertex;

	typedef _SPNode_iterator_global<DIM_TYPE, Cell, Data, Dim, Node&, Node*> iterator;
	typedef _SPNode_iterator_global<DIM_TYPE, Cell, Data, Dim, const Node&,
			const Node*> const_iterator;

	typedef _SPNodeFace_iterator_global<DIM_TYPE, Cell, Data, Dim, Face&, Face*> iterator_face;
	typedef _SPNodeFace_iterator_global<DIM_TYPE, Cell, Data, Dim, const Face&,
			const Face*> const_iterator_face;

	typedef _SPNodeVertex_iterator_global<DIM_TYPE, Cell, Data, Dim, Vertex&,
			Vertex*> iterator_vertex;
	typedef _SPNodeVertex_iterator_global<DIM_TYPE, Cell, Data, Dim,
			const Vertex&, const Vertex*> const_iterator_vertex;
protected:
	SpaceT<pTree, Dim> fm;
	SpaceT<int, Dim> att;

	void _delete();
public:

	//arrayListT<pNode>* pNode_Index;

	Forest() {
		//pNode_Index = NULL_PTR;
	}

	Forest(int i, int j, Float ox, Float oy, Float length, int maxl, int = 0,
			Float = 0);
	//Forest(int ilen, int jlen, Float ox, Float oy, Float length, int maxl,
	//		TreeFactory& tf);
	void reconstruct(size_type i, size_type j, size_type = 0);
	~Forest();

	void setpTree(pTree, int i, int j, int = 0);
	void set_attribution(int attr, int i, int j, int = 0);
	void set_attribution(int i, int en);
	int get_attribution(int i, int j, int = 0) const;
	int get_attribution_1d(int i) const {
		assert(i < att.size());
		return att.at_1d(i);
	}
//get pTree =================================
	size_type to_1d_idx(size_type i, size_type j, size_type k = 0) const {
		return fm.to_1d_idx(i, j, k);
	}

	pTree getpTree_1d(size_type i);
	const pTree getpTree_1d(size_type i) const;

	pTree getpTree(size_type i, size_type j, size_type = 0);
	const pTree getpTree(size_type i, size_type j, size_type = 0) const;

	pTree getpTree(const typename DIM_TYPE::Point&);
	pNode getpNode(const typename DIM_TYPE::Point&);
	//--------------------------
	pTree getMinXeTree() const;
	pTree getMaxXeTree() const;
	pTree getMinYeTree() const;
	pTree getMaxYeTree() const;
	pTree getMinZeTree() const;
	pTree getMaxZeTree() const;

	pTree getFirstpeTree() const;
	pTree getLastpeTree() const;

	size_type getFirstpeTreeIdx() const;
	size_type getLastpeTreeIdx() const;

	inline size_type iLen() const {
		return fm.iLen();
	}
	inline size_type jLen() const {
		return fm.jLen();
	}
	inline size_type kLen() const {
		return (Dim < 3) ? 0 : fm.kLen();
	}
	inline bool isEmpty() const {
		if (fm.size() <= 0) {
			return true;
		} else {
			return false;
		}
	}
	size_type get_dim() const {
		return Dim;
	}
	size_type size() const {
		return fm.size();
	}
	size_type getNumEnableTree() const {
		return att.count_equal(ATT_ENABLE);
	}
	size_type getNUmDisableTree() const {
		return att.count_equal(ATT_DISABLE);
	}
	Float getTreeDx() const;

//Iterator===========================
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;
	iterator last();
	const_iterator last() const;

	iterator_face begin_face();
	const_iterator_face begin_face() const;
	iterator_face end_face();
	const_iterator_face end_face() const;

	iterator_vertex begin_vertex();
	const_iterator_vertex begin_vertex() const;
	iterator_vertex end_vertex();
	const_iterator_vertex end_vertex() const;
//Index==============================
	//void IndexNodes();

//Connect==============================
	void ConnectTrees();

//count====================================
	int count_leaf() const;
	int check_leaf_with_data() const;
//IO=======================================
protected:
	void _draw_to_vtk_leaf(FILE*& file, int numleaf) const;
public:
	void show() const;
	void show_as_contour(LarusDef::size_type idx_) const;
	void show_info() const;
	//2D================================================
	void draw_to_gnuplot_root(std::string filename, int mode) const;
	void draw_to_gnuplot_leaf(std::string filename) const;
	//3D and 2D=========================================
	void draw_to_vtk_root(std::string filename) const;
	void draw_to_vtk_leaf(std::string filename) const;
//=========================================
};

//===============================================
typedef Forest<Dimension_3D> Forest3D;
typedef Forest<Dimension_2D> Forest2D;
//===============================================

//Function out side of class=====================
template<typename DIM_TYPE>
void toListpNode(Forest<DIM_TYPE>& forest,
		ListT<typename DIM_TYPE::pNode>& listp) {
	for (LarusDef::size_type i = 0; i < forest.iLen(); i++) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); j++) {
			for (LarusDef::size_type k = 0;
					k < ((DIM_TYPE::DIM == 3) ? forest.kLen() : 1); k++) {
				if (forest.get_attribution(i, j, k) == ATT_ENABLE)
					toList_leaf(forest.getpTree(i, j, k), listp);
			}
		}
	}
}

template<typename DIM_TYPE>
void toListpNode(Forest<DIM_TYPE>& forest, //forest
		ListT<typename DIM_TYPE::pNode>& listp, //list
		CSAxis aixs,  //aixs
		Float loc) {  //loc
	for (LarusDef::size_type i = 0; i < forest.iLen(); i++) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); j++) {
			for (LarusDef::size_type k = 0;
					k < ((DIM_TYPE::DIM == 3) ? forest.kLen() : 1); k++) {
				if (forest.getAttribution(i, j, k) == ATT_ENABLE)
					toList_leaf(forest.getpTree(i, j, k), listp, aixs, loc);
			}
		}
	}
}

template<typename DIM_TYPE>
void toListpNode(Forest<DIM_TYPE>& forest, //forest
		ListT<typename DIM_TYPE::pNode>& listp, //list
		CSAxis aixs,  //aixs
		Float loc, int level) {  //loc
	for (LarusDef::size_type i = 0; i < forest.iLen(); i++) {
		for (LarusDef::size_type j = 0; j < forest.jLen(); j++) {
			for (LarusDef::size_type k = 0;
					k < ((DIM_TYPE::DIM == 3) ? forest.kLen() : 1); k++) {
				if (forest.getAttribution(i, j, k) == ATT_ENABLE)
					toList_leaf(forest.getpTree(i, j, k), listp, aixs, loc,
							level);
			}
		}
	}
}

//===============================================

template<typename DIM_TYPE>
Forest<DIM_TYPE>::Forest(int ilen, int jlen, Float ox, Float oy, Float length,
		int maxl, int klen, Float oz) {
	//pNode_Index = NULL_PTR;
	if (Dim == 2) {  //
		fm.reconstruct(ilen, jlen);
		att.reconstruct(ilen, jlen);
		for (int i = 0; i < ilen; i++) {
			for (int j = 0; j < jlen; j++) {
				typename DIM_TYPE::Cell c(ox + length * i, oy + length * j,
						ox + (i + 1) * length, oy + (j + 1) * length, 0.0, 0.0);
				fm(i, j) = new Tree(c, maxl);
				att(i, j) = ATT_ENABLE;
			}
		}
	}
	if (Dim == 3) {
		fm.reconstruct(ilen, jlen, klen);
		att.reconstruct(ilen, jlen, klen);
		for (int i = 0; i < ilen; i++) {
			for (int j = 0; j < jlen; j++) {
				for (int k = 0; k < klen; k++) {
					typename DIM_TYPE::Cell c(
							//
							ox + length * i, oy + length * j, oz + length * k,
							ox + (i + 1) * length, oy + (j + 1) * length,
							oz + (k + 1) * length);
					fm(i, j, k) = new Tree(c, maxl);
					att(i, j, k) = ATT_ENABLE;
				}
			}
		}
	}
}
template<typename DIM_TYPE>
void Forest<DIM_TYPE>::_delete() {
	for (int i = 0; i < fm.size(); i++) {
		if (fm.at_1d(i) != NULL_PTR) {
			delete fm.at_1d(i);
		}
	}
}

template<typename DIM_TYPE>
Forest<DIM_TYPE>::~Forest() {
	//if (pNode_Index != NULL_PTR) {
	//delete pNode_Index;
	//}
	_delete();
}
template<typename DIM_TYPE>
void Forest<DIM_TYPE>::setpTree(Forest<DIM_TYPE>::pTree pt, int i, int j,
		int k) {
	if (Dim == 3 && i < fm.iLen() && i >= 0 && j < fm.jLen() && j >= 0
			&& k < fm.kLen() && k >= 0) {
		fm(i, j, k) = pt;
	} else if (Dim == 2 && i < fm.iLen() && i >= 0 && j < fm.jLen() && j >= 0) {
		fm(i, j) = pt;
	} else {
		std::cerr << "setpTree Error \n";
	}
}
template<typename DIM_TYPE>
void Forest<DIM_TYPE>::set_attribution(int i, int en) {
	if (i < att.size() && i >= 0) {
		att.at_1d(i) = en;
	}
}
template<typename DIM_TYPE>
void Forest<DIM_TYPE>::set_attribution(int attr, int i, int j, int k) {
	if (i < att.iLen() && i >= 0 && j < att.jLen() && j >= 0
			&& (Dim == 3 ? (k < att.kLen() && k >= 0) : true)) {
		att(i, j, k) = attr;
	}
}
template<typename DIM_TYPE>
int Forest<DIM_TYPE>::get_attribution(int i, int j, int k) const {
	if (i < att.iLen() && i >= 0 && j < att.jLen() && j >= 0
			&& (Dim == 3 ? (k < att.kLen() && k >= 0) : true)) {
		return att(i, j, k);
	}
	return ATT_DISABLE;
}
template<typename DIM_TYPE>
void Forest<DIM_TYPE>::reconstruct(Forest<DIM_TYPE>::size_type i,
		Forest<DIM_TYPE>::size_type j, Forest<DIM_TYPE>::size_type k) {
	fm.reconstruct(i, j, k);
	att.reconstruct(i, j, k);
	att.assign(ATT_ENABLE);
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::pTree Forest<DIM_TYPE>::getpTree_1d(size_type i) {
	if (i >= 0 && i < fm.size()) {
		if (att.at_1d(i) == ATT_ENABLE) {
			return fm.at_1d(i);
		} else {
			return NULL_PTR;
		}
	} else {
		return NULL_PTR;
	}
}

template<typename DIM_TYPE>
const typename Forest<DIM_TYPE>::pTree Forest<DIM_TYPE>::getpTree_1d(
		size_type i) const {
	if (i >= 0 && i < fm.size()) {
		if (att.at_1d(i) == ATT_ENABLE) {
			return fm.at_1d(i);
		} else {
			return NULL_PTR;
		}
	} else {
		return NULL_PTR;
	}
}

template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::pTree Forest<DIM_TYPE>::getpTree(
		Forest<DIM_TYPE>::size_type i, Forest<DIM_TYPE>::size_type j,
		Forest<DIM_TYPE>::size_type k) {
	if (fm.testIdxIJK(i, j, k)) {
		if (att(i, j, k) == ATT_ENABLE) {
			return fm(i, j, k);
		} else {
			return NULL_PTR;
		}
	} else {
		std::cerr << " >! Index out of range " << i << " " << j << " " << k
				<< "\n";
		return NULL_PTR;
	}
}
template<typename DIM_TYPE>
const typename Forest<DIM_TYPE>::pTree Forest<DIM_TYPE>::getpTree(
		Forest<DIM_TYPE>::size_type i, Forest<DIM_TYPE>::size_type j,
		Forest<DIM_TYPE>::size_type k) const {
	if (fm.testIdxIJK(i, j, k)) {
		if (att(i, j, k) == ATT_ENABLE) {
			return fm(i, j, k);
		} else {
			return NULL_PTR;
		}
	} else {
		std::cerr << " >! Index out of range " << i << " " << j << " " << k
				<< "\n";
		return NULL_PTR;
	}
}

template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::pTree Forest<DIM_TYPE>::getpTree(
		const typename DIM_TYPE::Point& p) {
	for (size_type i = 0; i < fm.iLen(); i++) {
		for (size_type j = 0; j < fm.jLen(); j++) {
			for (size_type k = 0; k < (Dim == 3 ? fm.kLen() : 1); k++) {
				if (att(i, j, k) == 1) {
					if (fm(i, j, k)->isInOnTree(p)) {
						return fm(i, j, k);
					}
				}
			}
		}
	}
	return NULL_PTR;
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::pNode Forest<DIM_TYPE>::getpNode(const typename DIM_TYPE::Point& p){
	pTree pt = getpTree(p);
	if(pt==NULL_PTR){
		return NULL_PTR;
	}else{
		return pt->getpNode(p);
	}
}

template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::pTree Forest<DIM_TYPE>::getMinXeTree() const {
	for (size_type i = 0; i < fm.iLen(); i++) {
		for (size_type j = 0; j < fm.jLen(); j++) {
			for (size_type k = 0; k < (Dim == 3 ? fm.kLen() : 1); k++) {
				if (att(i, j, k) == ATT_ENABLE)
					return fm(i, j, k);
			}
		}
	}
	return NULL_PTR;
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::pTree Forest<DIM_TYPE>::getMaxXeTree() const {
	for (size_type i = fm.iLen() - 1; i >= 0; i--) {
		for (size_type j = 0; j < fm.jLen(); j++) {
			for (size_type k = 0; k < (Dim == 3 ? fm.kLen() : 1); k++) {
				if (att(i, j, k) == ATT_ENABLE)
					return fm(i, j, k);
			}
		}
	}
	return NULL_PTR;
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::pTree Forest<DIM_TYPE>::getMinYeTree() const {
	for (size_type j = 0; j < fm.jLen(); j++) {
		for (size_type i = 0; i < fm.iLen(); i++) {
			for (size_type k = 0; k < (Dim == 3 ? fm.kLen() : 1); k++) {
				if (att(i, j, k) == ATT_ENABLE)
					return fm(i, j, k);
			}
		}
	}
	return NULL_PTR;
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::pTree Forest<DIM_TYPE>::getMaxYeTree() const {
	for (size_type j = fm.jLen() - 1; j >= 0; j--) {
		for (size_type i = 0; i < fm.iLen(); i++) {
			for (size_type k = 0; k < (Dim == 3 ? fm.kLen() : 1); k++) {
				if (att(i, j, k) == ATT_ENABLE)
					return fm(i, j, k);
			}
		}
	}
	return NULL_PTR;
}

template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::pTree Forest<DIM_TYPE>::getMinZeTree() const {
	ASSERT(Dim == 3);
	for (size_type k = 0; k < fm.kLen(); k++) {
		for (size_type i = 0; i < fm.iLen(); i++) {
			for (size_type j = 0; j < fm.jLen(); j++) {
				if (att(i, j, k) == ATT_ENABLE)
					return fm(i, j, k);
			}
		}
	}
	return NULL_PTR;
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::pTree Forest<DIM_TYPE>::getMaxZeTree() const {
	ASSERT(Dim == 3);
	for (size_type k = fm.kLen() - 1; k >= 0; k--) {
		for (size_type i = 0; i < fm.iLen(); i++) {
			for (size_type j = 0; j < fm.jLen(); j++) {
				if (att(i, j, k) == ATT_ENABLE)
					return fm(i, j, k);
			}
		}
	}
	return NULL_PTR;
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::pTree Forest<DIM_TYPE>::getFirstpeTree() const {
	for (size_type i = 0; i < fm.size(); ++i) {
		if (att.at_1d(i) == ATT_ENABLE) {
			return fm.at_1d(i);
		}
	}
	return NULL_PTR;
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::pTree Forest<DIM_TYPE>::getLastpeTree() const {
	for (size_type i = fm.size() - 1; i >= 0; --i) {
		if (att.at_1d(i) == ATT_ENABLE) {
			return fm.at_1d(i);
		}
	}
	return NULL_PTR;
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::size_type Forest<DIM_TYPE>::getFirstpeTreeIdx() const {
	for (size_type i = 0; i < fm.size(); ++i) {
		if (att.at_1d(i) == ATT_ENABLE) {
			return i;
		}
	}
	return -1;
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::size_type Forest<DIM_TYPE>::getLastpeTreeIdx() const {
	for (size_type i = fm.size() - 1; i >= 0; --i) {
		if (att.at_1d(i) == ATT_ENABLE) {
			return i;
		}
	}
	return -1;
}

template<typename DIM_TYPE>
Float Forest<DIM_TYPE>::getTreeDx() const {
	if (size() > 0) {
		pTree pt = NULL_PTR;
		for (size_type i = 0; i < size(); ++i) {
			if ((pt = this->getpTree_1d(i)) != NULL_PTR) {
				return pt->getpRootCell()->getDx();
			}
		}
	}
	return 0.0;
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::iterator Forest<DIM_TYPE>::begin() {
	size_type idx = this->getFirstpeTreeIdx();
	pTree pt = this->getFirstpeTree();
	pNode pn = getFirstLeaf(pt->getpRootNode());
	return iterator(this, idx, pn);
}

template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::const_iterator Forest<DIM_TYPE>::begin() const {
	size_type idx = this->getFirstpeTreeIdx();
	pTree pt = this->getFirstpeTree();
	pNode pn = getFirstLeaf(pt->getpRootNode());
	return iterator(this, idx, pn);
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::iterator Forest<DIM_TYPE>::last() {
	size_type idx = this->getLastpeTreeIdx();
	pTree pt = this->getLastpeTree();
	pNode pn = getLastLeaf(pt->getpRootNode());
	return iterator(this, idx, pn);
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::const_iterator Forest<DIM_TYPE>::last() const {
	size_type idx = this->getLastpeTreeIdx();
	pTree pt = this->getLastpeTree();
	pNode pn = getLastLeaf(pt->getpRootNode());
	return iterator(this, idx, pn);
}

template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::iterator Forest<DIM_TYPE>::end() {
	size_type idx = this->getLastpeTreeIdx();
	pTree pt = this->getLastpeTree();
	pNode pn = pt->getpRootNode();
	return iterator(this, idx, pn);
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::const_iterator Forest<DIM_TYPE>::end() const {
	size_type idx = this->getLastpeTreeIdx();
	pTree pt = this->getLastpeTree();
	pNode pn = pt->getpRootNode();
	return iterator(this, idx, pn);
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::iterator_face Forest<DIM_TYPE>::begin_face() {
	size_type idx = this->getFirstpeTreeIdx();
	pTree pt = this->getFirstpeTree();
	pNode pn = getFirstLeaf(pt->getpRootNode());
	pNode pnei = pn->getNeighborFast(SPD_IM);
	Face f(pn, pnei, SPD_IM, getFaceType(pn, pnei));
	return iterator_face(this, idx, f);
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::const_iterator_face Forest<DIM_TYPE>::begin_face() const {
	size_type idx = this->getFirstpeTreeIdx();
	pTree pt = this->getFirstpeTree();
	pNode pn = getFirstLeaf(pt->getpRootNode());
	pNode pnei = pn->getNeighborFast(SPD_IM);
	Face f(pn, pnei, SPD_IM, getFaceType(pn, pnei));
	return iterator_face(this, idx, f);
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::iterator_face Forest<DIM_TYPE>::end_face() {
	size_type idx = this->getLastpeTreeIdx();
	pTree pt = this->getLastpeTree();
	pNode pn = pt->getpRootNode();
	Face f(pn, NULL_PTR, SPD_IM, SPFT_Boundary);
	return iterator_face(this, idx, f);
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::const_iterator_face Forest<DIM_TYPE>::end_face() const {
	size_type idx = this->getLastpeTreeIdx();
	pTree pt = this->getLastpeTree();
	pNode pn = pt->getpRootNode();
	Face f(pn, NULL_PTR, SPD_IM, SPFT_Boundary);
	return iterator_face(this, idx, f);
}

template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::iterator_vertex Forest<DIM_TYPE>::begin_vertex() {
	size_type idx = this->getFirstpeTreeIdx();
	pTree pt = this->getFirstpeTree();
	pNode pn = getFirstLeaf(pt->getpRootNode());
	Vertex v(pn, SPD_MP);
	return iterator_vertex(this, idx, v);
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::const_iterator_vertex Forest<DIM_TYPE>::begin_vertex() const {
	size_type idx = this->getFirstpeTreeIdx();
	pTree pt = this->getFirstpeTree();
	pNode pn = getFirstLeaf(pt->getpRootNode());
	Vertex v(pn, SPD_MP);
	return const_iterator_vertex(this, idx, v);
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::iterator_vertex Forest<DIM_TYPE>::end_vertex() {
	size_type idx = this->getLastpeTreeIdx();
	pTree pt = this->getLastpeTree();
	pNode pn = pt->getpRootNode();
	Vertex v(pn, SPD_MP);
	return iterator_vertex(this, idx, v);
}
template<typename DIM_TYPE>
typename Forest<DIM_TYPE>::const_iterator_vertex Forest<DIM_TYPE>::end_vertex() const {
	size_type idx = this->getLastpeTreeIdx();
	pTree pt = this->getLastpeTree();
	pNode pn = pt->getpRootNode();
	Vertex v(pn, SPD_MP);
	return const_iterator_vertex(this, idx, v);
}

template<typename DIM_TYPE>
void Forest<DIM_TYPE>::ConnectTrees() {
	for (int i = 0; i < fm.iLen(); i++) {
		for (int j = 0; j < fm.jLen(); j++) {
			for (int k = 0; k < (Dim == 3 ? fm.kLen() : 1); k++) {
				if (att(i, j, k) == 1) {
					typename DIM_TYPE::pTree xp = NULL_PTR, xm = NULL_PTR, //x
							yp = NULL_PTR, ym = NULL_PTR, //y
							zp = NULL_PTR, zm = NULL_PTR; //z
					//E
					if (att.testIdx(0, i + 1)) {
						if (att(i + 1, j, k) == 1) {
							xp = fm(i + 1, j, k);
						}
					}
					//W
					if (att.testIdx(0, i - 1)) {
						if (att(i - 1, j, k) == 1) {
							xm = fm(i - 1, j, k);
						}
					}
					if (att.testIdx(1, j + 1)) {
						if (att(i, j + 1, k) == 1) {
							yp = fm(i, j + 1, k);
						}
					}
					if (att.testIdx(1, j - 1)) {
						if (att(i, j - 1, k) == 1) {
							ym = fm(i, j - 1, k);
						}
					}
					if (Dim == 3) {
						if (att.testIdx(2, k + 1)) {
							if (att(i, j, k + 1) == 1) {
								zp = fm(i, j, k + 1);
							}
						}
						if (att.testIdx(2, k - 1)) {
							if (att(i, j, k - 1) == 1) {
								zm = fm(i, j, k - 1);
							}
						}
					}
					fm(i, j, k)->setNeighborTree(xm, xp, ym, yp, zm, zp);
				}
			}
		}
	}
	for (int i = 0; i < fm.iLen(); i++) {
		for (int j = 0; j < fm.jLen(); j++) {
			for (int k = 0; k < (Dim == 3 ? fm.kLen() : 1); k++) {
				if (att(i, j, k) == 1) {
					fm(i, j, k)->ConnectNodes();
				}
			}
		}
	}
}
template<typename DIM_TYPE>
int Forest<DIM_TYPE>::count_leaf() const {
	int res = 0;
	if (!isEmpty()) {
		for (size_type i = 0; i < fm.size(); ++i) {
			if (att.at_1d(i) == ATT_ENABLE) {
				res += fm.at_1d(i)->count_leaf();
			}
		}
	}
	return res;
}

void is_leaf_has_data(pQTNode pn, utPointer p);
void is_leaf_has_data(pOCNode pn, utPointer p);

template<typename DIM_TYPE>
int Forest<DIM_TYPE>::check_leaf_with_data() const {
	int res = 0;
	if (!isEmpty()) {
		for (int i = 0; i < fm.iLen(); i++) {
			for (int j = 0; j < fm.jLen(); j++) {
				for (int k = 0; k < (Dim == 3 ? fm.kLen() : 1); k++) {
					if (att(i, j, k) == 1) {
						fm(i, j, k)->Traversal(is_leaf_has_data, &res);
					}
				}
			}
		}
	}
	int leaf = count_leaf();
	if (res == leaf) {
		std::cout << " > All leaves has data!  Check OK! \n";
	} else {
		std::cout << " > " << leaf - res << " leaves don't have data! \n";
	}
	return res;
}

template<typename DIM_TYPE>
void Forest<DIM_TYPE>::show() const {
	ASSERT(DIM_TYPE::DIM == 2);
	ListT<Segment2D> lseg;
	typedef typename DIM_TYPE::size_type st;
	for (st i = 0; i < fm.size(); ++i) {
		if (att.at_1d(i) == ATT_ENABLE) {
			typename DIM_TYPE::pCell cell = fm.at_1d(i)->getpRootNode()->cell;
			lseg.push_back(
					Segment2D(cell->getPoint(eCPL_M, eCPL_M),
							cell->getPoint(eCPL_M, eCPL_P)));
			lseg.push_back(
					Segment2D(cell->getPoint(eCPL_M, eCPL_P),
							cell->getPoint(eCPL_P, eCPL_P)));
			lseg.push_back(
					Segment2D(cell->getPoint(eCPL_P, eCPL_P),
							cell->getPoint(eCPL_P, eCPL_M)));
			lseg.push_back(
					Segment2D(cell->getPoint(eCPL_P, eCPL_M),
							cell->getPoint(eCPL_M, eCPL_M)));
		}
	}
	//show on gnuplot
	Gnuplot gp("lines");
	arrayList arrx(lseg.size() * 2);
	arrayList arry(lseg.size() * 2);
	st i = 0;
	for (ListT<Segment2D>::iterator iter = lseg.begin(); iter != lseg.end();
			++iter) {
		arrx[i] = iter->PSX();
		arry[i] = iter->PSY();
		arrx[i + 1] = iter->PEX();
		arry[i + 1] = iter->PEY();
		i += 2;
	}
	Float maxx = arrx.findMax();
	Float maxy = arry.findMax();
	Float minx = arrx.findMin();
	Float miny = arry.findMin();
	gp.set_equal_ratio();
	for (st i = 0; i < fm.size(); ++i) {
		if (att.at_1d(i) == ATT_ENABLE) {
			typename DIM_TYPE::pCell cell = fm.at_1d(i)->getpRootNode()->cell;
			Point2D cp = cell->getCenterPoint();
			std::ostringstream ss;
			ss << "\"" << int(i) << " ( " << int(i) / int(jLen()) << ", "
					<< int(i) % int(jLen()) << " )\" at first " << cp.x
					<< ", first " << cp.y << " center";
			gp.set_label(ss.str());
		}
	}
	gp.set_xrange(minx - (maxx - minx) * 0.05, maxx + (maxx - minx) * 0.05);
	gp.set_yrange(miny - (maxy - miny) * 0.05, maxy + (maxy - miny) * 0.05);
	gp.set_xlabel("X");
	gp.set_ylabel("Y");
	gp.plot_2_jump(arrx, arry, 2);
}


template<typename DIM_TYPE>
void Forest<DIM_TYPE>::show_as_contour(LarusDef::size_type idx_) const {
	typedef typename DIM_TYPE::CellData::value_type vt;
	ListT<vt> lxc, lyc, lxm, lxp, lym, lyp, lval;
	for (typename Forest<DIM_TYPE>::const_iterator iter = this->begin();
			iter != this->end(); iter++) {
		const typename DIM_TYPE::Node* pnode = iter.get_pointer();
		lxc.push_back(pnode->cell->getCenterPoint().x);
		lyc.push_back(pnode->cell->getCenterPoint().y);
		lxm.push_back(pnode->cell->getMM().x);
		lxp.push_back(pnode->cell->getPP().x);
		lym.push_back(pnode->cell->getMM().y);
		lyp.push_back(pnode->cell->getPP().y);
		lval.push_back(pnode->data->aCenterData[idx_]);
	}
	typename ListT<vt>::const_iterator iter = lval.begin();
	vt max = (*iter);
	vt min = (*iter);
	for (++iter; iter != lval.end(); ++iter) {
		if ((*iter) > max) {
			max = (*iter);
		}
		if ((*iter) < min) {
			min = (*iter);
		}
	}
	Gnuplot gp("boxes");
	std::string cmdstr = "with boxxy title \"\" fs solid palette";
	gp.set_palette_blue_red();
	if (max == min) {
		gp.set_cbrange(max - 1, max);
	} else {
		gp.set_cbrange(min, max);
	}
	std::ostringstream ss;
	gp.set_xlabel(ss.str());
	gp.plot_7(lxc, lyc, lxm, lxp, lym, lyp, lval, cmdstr);
}

template<typename DIM_TYPE>
void Forest<DIM_TYPE>::show_info() const {
	std::cout << "=>Forest Info: <========\n";
	std::cout << "Dim          :" << Dim << "\n";
	std::cout << "Size         :" << fm.iLen() << " x " << fm.jLen();
	if (Dim == 3) {
		std::cout << " x " << fm.kLen() << "\n";
	} else {
		std::cout << "\n";
	}
	std::cout << "Num of Tree  :" << size() << "\n";
	std::cout << "Num of e Tree:" << getNumEnableTree() << "\n";
	//std::cout << "Num of leaf  :" << count_leaf() << std::endl;
	if (size() > 0) {
		std::cout << "Length       :" << getTreeDx() << std::endl;
		pTree pt = getMinXeTree();
		if (pt != NULL) {
			std::cout << "x e range min:"
					<< pt->getpRootCell()->getPoint(eCPL_M, eCPL_M, eCPL_M).x
					<< "\n";
		}
		pt = getMaxXeTree();
		if (pt != NULL) {
			std::cout << "          max:"
					<< pt->getpRootCell()->getPoint(eCPL_P, eCPL_P, eCPL_P).x
					<< "\n";
		}
		pt = getMinYeTree();
		if (pt != NULL) {
			std::cout << "y e range min:"
					<< pt->getpRootCell()->getPoint(eCPL_M, eCPL_M, eCPL_M).y
					<< "\n";
		}
		pt = getMaxYeTree();
		if (pt != NULL) {
			std::cout << "          max:"
					<< pt->getpRootCell()->getPoint(eCPL_P, eCPL_P, eCPL_P).y
					<< "\n";
		}
		if (Dim == 3) {
			pt = getMinZeTree();
			if (pt != NULL) {
				std::cout << "z e range min:"
						<< pt->getpRootCell()->getPoint(eCPL_M, eCPL_M, eCPL_M)[2]
						<< "\n";
			}
			pt = getMaxZeTree();
			if (pt != NULL) {
				std::cout << "          max:"
						<< pt->getpRootCell()->getPoint(eCPL_P, eCPL_P, eCPL_P)[2]
						<< "\n";
			}
		}
	}
	//Leaf information---------------------------
	_IF_TRUE_RETRUN(isEmpty());

	pTree pt = getMinXeTree();
	_IF_TRUE_RETRUN(pt==NULL_PTR);
	int _maxlevel = pt->getMaxLevel();
	arrayList_int num_node(_maxlevel + 1);
	arrayList_int num_leaf(_maxlevel + 1);
	for (int i = 0; i < size(); i++) {
		pTree tree = getpTree_1d(i);
		if (tree != NULL_PTR) {
			for (int il = 0; il <= _maxlevel; il++) {
				num_node[il] += tree->count_level(il);
				num_leaf[il] += tree->count_leaf_at_level(il);
			}
		}
	}
	std::cout << "level   Num   leaf   ratio%\n";
	int totalnode = 0;
	int totalleaf = 0;
	for (int i = 0; i <= _maxlevel; i++) {
		int nn = num_node[i];
		std::cout.flags(std::ios::right);
		std::cout.width(4);
		std::cout << i;
		std::cout.width(7);
		std::cout << nn;
		std::cout.width(7);
		int nl = num_leaf[i];
		std::cout << nl;
		std::cout.width(7);
		std::cout.precision(1);
		std::cout.setf(std::ios::fixed, std::ios::floatfield);
		std::cout << Float(nl) / Float(nn) * 100 << std::endl;
		totalnode += nn;
		totalleaf += nl;
	}
	std::cout.flags(std::ios::right);
	std::cout.width(11);
	std::cout << totalnode;
	std::cout.width(7);
	std::cout << totalleaf;
	std::cout.width(7);
	std::cout.precision(1);
	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout << Float(totalleaf) / Float(totalnode) * 100 << std::endl;
}
template<typename DIM_TYPE>
void Forest<DIM_TYPE>::draw_to_gnuplot_root(std::string filename,
		int mode) const {
	ASSERT(Dim == 2);
	FILE *data = open_file(filename, mode);
	for (int i = 0; i < fm.iLen(); i++) {
		for (int j = 0; j < fm.jLen(); j++) {
			if (this->get_attribution(i, j) == ATT_ENABLE) {
				typename DIM_TYPE::pCell cell =
						this->getpTree(i, j)->getpRootNode()->cell;
				fprintf(data, "%f %f %d %d %d\n",
						cell->getPoint(eCPL_M, eCPL_M).x,
						cell->getPoint(eCPL_M, eCPL_M).y, get_attribution(i, j),
						i, j);
				fprintf(data, "%f %f %d %d %d\n\n",
						cell->getPoint(eCPL_M, eCPL_P).x,
						cell->getPoint(eCPL_M, eCPL_P).y, get_attribution(i, j),
						i, j);

				fprintf(data, "%f %f %d %d %d\n",
						cell->getPoint(eCPL_M, eCPL_P).x,
						cell->getPoint(eCPL_M, eCPL_P).y, get_attribution(i, j),
						i, j);
				fprintf(data, "%f %f %d %d %d\n\n",
						cell->getPoint(eCPL_P, eCPL_P).x,
						cell->getPoint(eCPL_P, eCPL_P).y, get_attribution(i, j),
						i, j);

				fprintf(data, "%f %f %d %d %d\n",
						cell->getPoint(eCPL_P, eCPL_P).x,
						cell->getPoint(eCPL_P, eCPL_P).y, get_attribution(i, j),
						i, j);
				fprintf(data, "%f %f %d %d %d\n\n",
						cell->getPoint(eCPL_P, eCPL_M).x,
						cell->getPoint(eCPL_P, eCPL_M).y, get_attribution(i, j),
						i, j);

				fprintf(data, "%f %f %d %d %d\n",
						cell->getPoint(eCPL_P, eCPL_M).x,
						cell->getPoint(eCPL_P, eCPL_M).y, get_attribution(i, j),
						i, j);
				fprintf(data, "%f %f %d %d %d\n\n",
						cell->getPoint(eCPL_M, eCPL_M).x,
						cell->getPoint(eCPL_M, eCPL_M).y, get_attribution(i, j),
						i, j);
			}
		}
	}
}
template<typename DIM_TYPE>
void Forest<DIM_TYPE>::draw_to_gnuplot_leaf(std::string filename) const {
	ASSERT(Dim == 2);
	creat_empty_file(filename);
	for (size_type i = 0; i < fm.iLen(); i++) {
		for (size_type j = 0; j < fm.jLen(); j++) {
			const pTree pt = this->getpTree(i, j);
			if (pt != NULL_PTR)
				pt->draw_gnuplot(filename, 2);
		}
	}
}

template<typename DIM_TYPE>
void Forest<DIM_TYPE>::draw_to_vtk_root(std::string filename) const {
	//arrayList arr;
	size_type numleaf = this->getNumEnableTree();
	size_type vertexes = DIM_TYPE::Cell::NUM_VERTEXES;
	FILE *data = open_file(filename, 1);
	vtk_unstructured_grid_head(data, "Forest roots");
	fprintf(data, "POINTS %d float\n", numleaf * vertexes);
	for (size_type i = 0; i < fm.iLen(); i++) {
		for (size_type j = 0; j < fm.jLen(); j++) {
			for (size_type k = 0; k < (Dim == 3 ? fm.kLen() : 1); k++) {
				if (get_attribution(i, j, k) == ATT_ENABLE)
					fm(i, j, k)->getpRootCell()->output_vertex_in_vtk_order(
							data);
			}
		}
	}

	fprintf(data, "\n");
	fprintf(data, "CELLS %d %d \n", numleaf, numleaf * (vertexes + 1));

	for (int i = 0; i < numleaf * vertexes; i++) {
		if ((i % vertexes) == 0) {
			fprintf(data, "%d ", vertexes);
		}
		fprintf(data, "%d ", i);
		if ((i % vertexes) == vertexes - 1) {
			fprintf(data, "\n");
		}
	}
	fprintf(data, "\n");
	fprintf(data, "CELL_TYPES %d\n", numleaf);
	for (int i = 0; i < numleaf; i++) {
		if (Dim == 3) {
			fprintf(data, "%d \n", 11);
		} else if (Dim == 2) {
			fprintf(data, "%d \n", 8);
		}
	}
	fclose(data);
}

template<typename DIM_TYPE>
void Forest<DIM_TYPE>::draw_to_vtk_leaf(std::string filename) const {
	//arrayList arr;
	size_type numleaf = count_leaf();
	FILE *data = open_file(filename, 1);
	vtk_unstructured_grid_head(data, "Forest leafs");
	_draw_to_vtk_leaf(data, numleaf);
	fclose(data);
}

template<typename DIM_TYPE>
void Forest<DIM_TYPE>::_draw_to_vtk_leaf(FILE*& data, int numleaf) const {
	size_type vertexes = DIM_TYPE::Cell::NUM_VERTEXES;
	fprintf(data, "POINTS %d float\n", numleaf * vertexes);
	for (size_type i = 0; i < fm.iLen(); i++) {
		for (size_type j = 0; j < fm.jLen(); j++) {
			for (size_type k = 0; k < (Dim == 3 ? fm.kLen() : 1); k++) {
				if (get_attribution(i, j, k) == ATT_ENABLE)
					fm(i, j, k)->draw_to_vtk(data);
			}
		}
	}
	fprintf(data, "\n");
	fprintf(data, "CELLS %d %d \n", numleaf, numleaf * (vertexes + 1));

	for (int i = 0; i < numleaf * vertexes; i++) {
		if ((i % vertexes) == 0) {
			fprintf(data, "%d ", vertexes);
		}
		fprintf(data, "%d ", i);
		if ((i % vertexes) == vertexes - 1) {
			fprintf(data, "\n");
		}
	}
	fprintf(data, "\n");
	fprintf(data, "CELL_TYPES %d\n", numleaf);
	for (int i = 0; i < numleaf; i++) {
		if (Dim == 3) {
			fprintf(data, "%d \n", 11);
		} else if (Dim == 2) {
			fprintf(data, "%d ", 8);
		}
	}
	fprintf(data, "\n");
}

//===============================================
template<class Dimension>
SPTree<typename Dimension::Node, Dimension::DIM>* getNextpeTree(
		const Forest<Dimension>* f, LarusDef::size_type& i) {
	if (i >= 0 && i < f->size()) {
		LarusDef::size_type count = 0;
		for (LarusDef::size_type ii = (i + 1) % f->size(); count < f->size();
				++ii, ++count) {
			SPTree<typename Dimension::Node, Dimension::DIM>* pt =
					f->getpTree_1d(ii);
			if (NULL_PTR != pt) {
				i = ii;
				return pt;

			}
		}
		return NULL_PTR;
	} else {
		return NULL_PTR;
	}
}

template<class Dimension>
SPTree<typename Dimension::Node, Dimension::DIM>* getPrevpeTree(
		const Forest<Dimension>* f, LarusDef::size_type& i) {
	if (i >= 0 && i < f->getSize()) {
		LarusDef::size_type count = 0;
		LarusDef::size_type len = f->getSize();
		for (LarusDef::size_type ii = ((i - 1) < 0 ? (len - 1) : (i - 1))
				% f->getSize(); count < f->getSize(); --ii, ++count) {
			SPTree<typename Dimension::Node, Dimension::DIM>* pt =
					f->getpTree_1d(ii);
			if (NULL_PTR != pt) {
				i = ii;
				return pt;
			}
		}
		return NULL_PTR;
	} else {
		return NULL_PTR;
	}
}

//SPNode_iterator====================================
template<class Dimension, class Cell, class Data, LarusDef::size_type Dim,
		class _Ref, class _Ptr>
class _SPNode_iterator_global {
public:
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;
	typedef forward_iterator_tag iterator_category;

	typedef SPNode<Cell, Data, Dim> _Node;
	typedef Forest<Dimension> _Forest;
	typedef SPTree<_Node, Dim> _Tree;

	typedef _SPNode_iterator_global<Dimension, Cell, Data, Dim, _Node&, _Node*> iterator;
	typedef _SPNode_iterator_global<Dimension, Cell, Data, Dim, const _Node&,
			const _Node*> const_iterator;
	typedef _SPNode_iterator_global<Dimension, Cell, Data, Dim, _Ref, _Ptr> _Self;

	typedef _Node value_type;
	typedef _Ptr pointer;
	typedef _Ref reference;

	const _Forest* _f;
	LarusDef::size_type _idx;
	_Node* _ptr;

	_SPNode_iterator_global() {
		_ptr = NULL_PTR;
		_f = NULL_PTR;
		_idx = 0;
	}
	_SPNode_iterator_global(const _Forest* _f, LarusDef::size_type _idx,
			_Node* _ptr) {
		this->_ptr = _ptr;
		this->_f = _f;
		this->_idx = _idx;
	}
	_SPNode_iterator_global(const iterator& _x) {
		this->_ptr = _x._ptr;
		this->_f = _x._f;
		this->_idx = _x._idx;
	}

	void _incr() {
		_Node* end = _f->getLastpeTree()->getpRootNode();
		if (_ptr == end) {
			return;  //will not increase
		}
		_Node* s = getSiblingPlus(_ptr);
		if (s->father != NULL_PTR) {
			_ptr = getFirstLeaf(s);
		} else {
			if (s == end) {
				_ptr = s;
			} else {
				_Tree* pt = getNextpeTree(_f, _idx); //this will change _idx
				_ptr = getFirstLeaf(pt->getpRootNode());
			}
		}
	}

	bool operator==(const _SPNode_iterator_global& _x) const {
		return _ptr == _x._ptr;
	}
	bool operator!=(const _SPNode_iterator_global& _x) const {
		return _ptr != _x._ptr;
	}

	reference operator*() const {
		return (*_ptr);
	}

	pointer operator->() const {
		return &(operator*());
	}

	_Self & operator++() {
		this->_incr();
		return *this;
	}

	_Self operator++(int) {
		_Self __tmp = *this;
		this->_incr();
		return __tmp;
	}

	bool isExist() {
		return _ptr != NULL_PTR;
	}

	pointer get_pointer() {
		return _ptr;
	}

	const pointer get_pointer() const{
		return _ptr;
	}
};

//SPFace iterator ===============================
template<class Dimension, class Cell, class Data, LarusDef::size_type Dim,
		class _Ref, class _Ptr>
class _SPNodeFace_iterator_global {
public:
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;
	typedef forward_iterator_tag iterator_category;

	typedef SPNode<Cell, Data, Dim> _Node;
	typedef SPNodeFace<Cell, Data, Dim> _Face;
	typedef Forest<Dimension> _Forest;
	typedef SPTree<_Node, Dim> _Tree;

	typedef _SPNodeFace_iterator_global<Dimension, Cell, Data, Dim, _Face&,
			_Face*> iterator;
	typedef _SPNodeFace_iterator_global<Dimension, Cell, Data, Dim,
			const _Face&, const _Face*> const_iterator;
	typedef _SPNodeFace_iterator_global<Dimension, Cell, Data, Dim, _Ref, _Ptr> _Self;

	typedef _Node value_type;
	typedef _Ptr pointer;
	typedef _Ref reference;

	const _Forest* _f;
	LarusDef::size_type _idx;
	_Face _face;

	_SPNodeFace_iterator_global() {
		_face = _Face();
		_f = NULL_PTR;
		_idx = 0;
	}
	_SPNodeFace_iterator_global(const _Forest* _f, LarusDef::size_type _idx,
			const _Face& _pfa) {
		this->_face = _pfa;
		this->_f = _f;
		this->_idx = _idx;
	}
	_SPNodeFace_iterator_global(const iterator& _x) {
		this->_face = _x._face;
		this->_f = _x._f;
		this->_idx = _x._idx;
	}

	void _incr() {
		_Node* end = _f->getLastpeTree()->getpRootNode();
		if (_face.pnode == end) {
			return;  //will not increase
		}
		if (_face.direction < ((Dim == 2) ? 7 : 9)) {
			_face.direction = toDirection(int(_face.direction) + 1);
			_face.pneighbor = _face.pnode->getNeighborFast(_face.direction);
			_face.face_type = getFaceType(_face.pnode, _face.pneighbor);
			return;
		} else {
			_Node* s = getSiblingPlus(_face.pnode);
			if (s->father != NULL_PTR) {
				_face.pnode = getFirstLeaf(s);
			} else {
				if (s == end) {
					_face.pnode = s;
				} else {
					_Tree* pt = getNextpeTree(_f, _idx); //this will change _idx
					_face.pnode = getFirstLeaf(pt->getpRootNode());
				}
			}
			if (_face.pnode != end) {
				_face.direction = SPD_IM;
				_face.pneighbor = _face.pnode->getNeighborFast(SPD_IM);
				_face.face_type = getFaceType(_face.pnode, _face.pneighbor);
			} else {            //the last face-----------------
				_face.direction = SPD_IM;
				_face.pneighbor = NULL;
				_face.face_type = SPFT_Boundary;
			}
		}
	}

	bool operator==(const _SPNodeFace_iterator_global& _x) const { // e
		return _face == _x._face;
	}
	bool operator!=(const _SPNodeFace_iterator_global& _x) const { //e
		return _face != _x._face;
	}

	reference operator*() {
		return (_face);
	}

	pointer operator->() {
		return &(_face);
	}

	_Self & operator++() {
		this->_incr();
		return *this;
	}

	_Self operator++(int) {
		_Self __tmp = *this;
		this->_incr();
		return __tmp;
	}

	bool isExist() {
		return _face != NULL_PTR;
	}

	pointer get_pointer() {
		return &_face;
	}

	LarusDef::size_type get_TreeIdx() {
		return _idx;
	}

	const _Forest* get_pForest() {
		return _f;
	}
};

//SPNodeVertex iterator ===============================
template<class Dimension, class Cell, class Data, LarusDef::size_type Dim,
		class _Ref, class _Ptr>
class _SPNodeVertex_iterator_global {
public:
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;
	typedef forward_iterator_tag iterator_category;

	typedef SPNode<Cell, Data, Dim> _Node;
	typedef SPNodeVertex<Cell, Data, Dim> _Vertex;
	typedef Forest<Dimension> _Forest;
	typedef SPTree<_Node, Dim> _Tree;

	typedef _SPNodeVertex_iterator_global<Dimension, Cell, Data, Dim, _Vertex&,
			_Vertex*> iterator;
	typedef _SPNodeVertex_iterator_global<Dimension, Cell, Data, Dim,
			const _Vertex&, const _Vertex*> const_iterator;
	typedef _SPNodeVertex_iterator_global<Dimension, Cell, Data, Dim, _Ref, _Ptr> _Self;

	typedef _Node value_type;
	typedef _Ptr pointer;
	typedef _Ref reference;

	const _Forest* _f;
	LarusDef::size_type _idx;
	_Vertex _vertex;

	_SPNodeVertex_iterator_global() {
		_vertex = _Vertex();
		_f = NULL_PTR;
		_idx = 0;
	}
	_SPNodeVertex_iterator_global(const _Forest* _f, LarusDef::size_type _idx,
			const _Vertex& _pfa) {
		this->_vertex = _pfa;
		this->_f = _f;
		this->_idx = _idx;
	}
	_SPNodeVertex_iterator_global(const iterator& _x) {
		this->_vertex = _x._vertex;
		this->_f = _x._f;
		this->_idx = _x._idx;
	}

	void _incr() {
		_Node* end = _f->getLastpeTree()->getpRootNode();
		if (_vertex.pnode == end) {
			return;  //will not increase
		}
		if (_vertex.direction < ((Dim == 2) ? 3 : 25)) {
			_vertex.direction = toDirection(int(_vertex.direction) + 1);
			return;
		} else {
			_Node* s = getSiblingPlus(_vertex.pnode);
			if (s->father != NULL_PTR) {
				_vertex.pnode = getFirstLeaf(s);
			} else {
				if (s == end) {
					_vertex.pnode = s;
				} else {
					_Tree* pt = getNextpeTree(_f, _idx); //this will change _idx
					_vertex.pnode = getFirstLeaf(pt->getpRootNode());
				}
			}
			if (_vertex.pnode != end) {
				_vertex.direction = ((Dim == 2) ? SPD_MP : SPD_MMM);
			} else { //the last face-----------------
				_vertex.direction = ((Dim == 2) ? SPD_MP : SPD_MMM);
			}
		}
	}

	bool operator==(const _SPNodeVertex_iterator_global& _x) const { // e
		return _vertex == _x._vertex;
	}
	bool operator!=(const _SPNodeVertex_iterator_global& _x) const { //e
		return _vertex != _x._vertex;
	}

	reference operator*() {
		return (_vertex);
	}

	pointer operator->() {
		return &(_vertex);
	}

	_Self& operator++() {
		this->_incr();
		return (*this);
	}

	_Self operator++(int) {
		_Self __tmp = *this;
		this->_incr();
		return __tmp;
	}

	bool isExist() {
		return _vertex != NULL_PTR;
	}

	pointer get_pointer() {
		return _vertex;
	}
};

}

#endif /* FOREST_H_ */
