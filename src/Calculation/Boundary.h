/*
 * Boundary.h
 *
 *  Created on: May 9, 2015
 *      Author: zhou
 */

#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "../TypeDef.h"
#include "CalDef.h"
#include "../Geometry/Segment.h"
#include "../Grid/SPTree.h"
#include "../Grid/SPTreeNode.h"
#include "../Grid/Forest.h"
#include "../Algebra/Arithmetic.h"
#include "../Algebra/Expression.h"
#include "../IO/Gnuplot.h"
#include <map>
#include <set>

namespace Larus {

enum BoundaryType {
	BT_DIRICHLET = 1, BT_NEUMANN = 2, BT_SYMMETRIC = 3
};

inline bool isBoundary(pQuadTree pt, SPDirection dir) {
	ASSERT(pt != NULL_PTR);
	return NULL_PTR == pt->getNeighborpTree(dir);
}

inline bool isBoundary(pOCTree pt, SPDirection dir) {
	ASSERT(pt != NULL_PTR);
	return NULL_PTR == pt->getNeighborpTree(dir);
}

inline bool isBoundary(pQTNode pt, SPDirection dir) {
	ASSERT(pt != NULL_PTR);
	return NULL_PTR == pt->getNeighborFast(dir);
}

inline bool isBoundary(pOCNode pt, SPDirection dir) {
	ASSERT(pt != NULL_PTR);
	return NULL_PTR == pt->getNeighborFast(dir);
}

inline bool isBoundaryNode(pQTNode pt) {
	ASSERT(pt != NULL_PTR);
	const SPDirection FD[] = { SPD_IM, SPD_JP, SPD_IP, SPD_JM };
	for (int i = 0; i < 4; ++i) {
		if (isBoundary(pt, FD[i])) {
			return true;
		}
	}
	return false;
}

inline bool isBoundaryNode(pOCTree pt, SPDirection dir) {
	ASSERT(pt != NULL_PTR);
	return NULL_PTR == pt->getNeighborpTree(dir);
}

template<class DIMENSION>
struct BoundaryCondition {
	typedef int (*pFun_Boundary_set)(const BoundaryCondition&,
			typename DIMENSION::pNode //pGostNode
			);
	int tree_idx;
	SPDirection direction;
	LarusDef::size_type value_idx;
	utPointer utp;
	pFun_Boundary_set pfun;

	void show() const{
		std::cout<<"BC -----\n";
		std::cout<<"tree id  :"<<tree_idx<<"\n";
		std::cout<<"direction:"<<int(direction)<<"\n";
		std::cout<<"value id :"<<value_idx<<"\n";
	}
};

template<class DIMENSION>
struct BC_compare {
	bool operator()(const BoundaryCondition<DIMENSION>& lhs,
			const BoundaryCondition<DIMENSION>& rhs) const {
		if (lhs.tree_idx < rhs.tree_idx) {
			return true;
		} else if (lhs.tree_idx == rhs.tree_idx) {
			if (int(lhs.direction) < int(rhs.direction)) {
				return true;
			} else if (int(lhs.direction) == int(rhs.direction)) {
				return lhs.value_idx < rhs.value_idx;
			}
		}
		return false;
	}
};

template<class DIMENSION>
struct GhostID {
	int tree_idx;
	int node_idx;
	typename DIMENSION::pNode pnode;
	SPDirection direction;
};

template<class DIMENSION>
struct GhostID_compare {
	bool operator()(const GhostID<DIMENSION>& lhs,
			const GhostID<DIMENSION>& rhs) const {
		if (lhs.node_idx < rhs.node_idx) {
			return true;
		} else if (lhs.node_idx == rhs.node_idx) {
			return int(lhs.direction) < int(rhs.direction);
		} else {
			return false;
		}
	}
};

template<class DIMENSION>
class BCManager: public ObjectBase {
public:
	typedef typename DIMENSION::size_type st;
	typedef typename DIMENSION::pTree pTree;
	typedef typename DIMENSION::Tree Tree;
	typedef typename DIMENSION::pNode pNode;
	typedef typename DIMENSION::Node Node;

	typedef std::pair<GhostID<DIMENSION>, pNode> GhostNode;
	typedef std::map<GhostID<DIMENSION>, pNode, GhostID_compare<DIMENSION> > GhostMap;
	typedef std::set<BoundaryCondition<DIMENSION>, BC_compare<DIMENSION> > BCSet;
protected:
	BCSet bcset;
	GhostMap ghostmap;
public:
	Forest<DIMENSION>* pforest;

	BCManager(Forest<DIMENSION>* pf) :
			bcset() {
		pforest = pf;
	}

	~BCManager() {
		delete_ghost_nodes();
	}

	//new and delete  ===========================
	void new_ghost_nodes();
	void delete_ghost_nodes();

	//add bc condition ==========================
	void add_BC(const BoundaryCondition<DIMENSION>& bc);

	//find ghost ================================
	pNode find_ghost(int nid, SPDirection dir);

	//set boundary condition ====================
	void set_bc_on_ghost_nodes( //
			typename GhostMap::iterator&, typename BCSet::iterator&);
	void set_bc();
	//output=====================================
	void show_info() const;
	void show_set_of_boundary_condition() const;  //
	void show_tree_boundary() const;
	void show_boundary_contour(st v_idx) const;

	void draw_boundary_node(const std::string& filename) const;
	void draw_ghost_node(const std::string& filename) const;
};

template<class DIMENSION>
void BCManager<DIMENSION>::new_ghost_nodes() {
	for (Forest2D::iterator_face iterf = pforest->begin_face();
			iterf != pforest->end_face(); ++iterf) {
		if (iterf->face_type == SPFT_Boundary) {
			GhostID<DIMENSION> gid;
			gid.tree_idx = iterf.get_TreeIdx();
			gid.pnode = iterf->pnode;
			gid.node_idx = iterf->pnode->data->aCenterData[Idx_IDX];
			gid.direction = iterf->direction;
			pNode pn = new_ghost_node((*iterf));
			//change the index id ------------------------
			pn->data->aCenterData[Idx_IDX] = -gid.node_idx -1; //negative
			ghostmap.insert(GhostNode(gid, pn));
		}
	}
}
template<class DIMENSION>
void BCManager<DIMENSION>::delete_ghost_nodes() {
	for (typename GhostMap::iterator it = ghostmap.begin();
			it != ghostmap.end(); ++it) {
		if (it->second != NULL_PTR) {
			delete it->second;
			it->second = NULL_PTR;
		}
	}
}

template<class DIMENSION>
void BCManager<DIMENSION>::show_boundary_contour(
		typename BCManager<DIMENSION>::st v_idx) const {
	ASSERT(DIMENSION::DIM == 2);
	typedef typename DIMENSION::CellData::value_type vt;
	ListT<vt> lxc, lyc, lxm, lxp, lym, lyp, lval;
	for (typename GhostMap::const_iterator iter = ghostmap.begin();
			iter != ghostmap.end(); ++iter) {
		pNode pnode = iter->second;
		lxc.push_back(pnode->cell->getCenterPoint().x);
		lyc.push_back(pnode->cell->getCenterPoint().y);
		lxm.push_back(pnode->cell->getMM().x);
		lxp.push_back(pnode->cell->getPP().x);
		lym.push_back(pnode->cell->getMM().y);
		lyp.push_back(pnode->cell->getPP().y);
		lval.push_back(pnode->data->aCenterData[v_idx]);
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
		gp.set_cbrange(ABS(max), -ABS(max));
	} else {
		gp.set_cbrange(min, max);
	}

	std::ostringstream ss;
	gp.set_xlabel(ss.str());
	gp.set_equal_ratio();
	gp.plot_7(lxc, lyc, lxm, lxp, lym, lyp, lval, cmdstr);

}

template<class DIMENSION>
void BCManager<DIMENSION>::show_tree_boundary() const {
	ASSERT(DIMENSION::DIM == 2);
	ListT<Segment2D> lseg;
	pTree pt = NULL_PTR;
	for (st i = 0; i < pforest->size(); ++i) {
		pt = pforest->getpTree_1d(i);
		if (pt != NULL_PTR) {
			// check neighbor on all face direction
			// direction 4 5 6 7
			for (st i = 4; i <= 7; ++i) {
				pTree nt = pt->getNeighborpTree(toDirection(i));
				SPDirection d1, d2;
				if (nt == NULL_PTR) {
					if (i == 4 || i == 6) {
						d1 = Direction_Compose(toDirection(i), SPD_JM);
						d2 = Direction_Compose(toDirection(i), SPD_JP);
					} else {
						d1 = Direction_Compose(SPD_IM, toDirection(i));
						d2 = Direction_Compose(SPD_IP, toDirection(i));
					}
					Point2D p1 = pt->getpRootNode()->getPoint(d1);
					Point2D p2 = pt->getpRootNode()->getPoint(d2);
					if (i == 6 || i == 7) { //keep segment direction;
						Segment2D seg;
						seg.reconstruct(p1, p2);
						lseg.push_back(seg);
					} else {
						Segment2D seg;
						seg.reconstruct(p2, p1);
						lseg.push_back(seg);
					}
				}
			}
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
	gp.set_xrange(minx - (maxx - minx) * 0.05, maxx + (maxx - minx) * 0.05);
	gp.set_yrange(miny - (maxy - miny) * 0.05, maxy + (maxy - miny) * 0.05);
	gp.set_xlabel("X");
	gp.set_ylabel("Y");
	gp.plot_2_jump(arrx, arry, 2);
}

template<class DIMENSION>
void BCManager<DIMENSION>::draw_boundary_node(const std::string& filename) const {
	ASSERT(DIMENSION::DIM == 2);
	FILE* file = open_file(filename, 1);
	for (Forest2D::iterator iter = pforest->begin(); iter != pforest->end();
			++iter) {
		if (isBoundaryNode(iter.get_pointer())) {
			draw_gnuplot_boundary(iter.get_pointer(), file, SP_XY);
		}
	}
	fclose(file);
}

template<class DIMENSION>
void BCManager<DIMENSION>::draw_ghost_node(const std::string& filename) const {
	ASSERT(DIMENSION::DIM == 2);
	FILE* file = open_file(filename, 1);
	for (typename GhostMap::const_iterator it = ghostmap.begin();
			it != ghostmap.end(); ++it) {
		draw_gnuplot_boundary(it->second, file, SP_XY);
	}
	fclose(file);
}
template<class DIMENSION>
void BCManager<DIMENSION>::show_set_of_boundary_condition() const {
	std::cout << "BC set:  size"<<bcset.size()<<"\n";
	std::cout << "tree idx   direction   value idx \n";
	for (typename BCSet::const_iterator iter = bcset.begin();
			iter != bcset.end(); ++iter) {
		std::cout << iter->tree_idx << "    " << iter->direction << "    "
				 << iter->value_idx << " \n";
	}
}

template<class DIMENSION>
void BCManager<DIMENSION>::add_BC(const BoundaryCondition<DIMENSION>& bc) {
	std::pair<typename BCSet::iterator,bool> ret;
	ret = bcset.insert(bc);
	if( ret.second == false ){
		std::cout<<"Insert BC failed! ---------";
		bc.show();
	}
}

template<class DIMENSION>
typename BCManager<DIMENSION>::pNode BCManager<DIMENSION>::find_ghost(int nid,
		SPDirection dir) {
	GhostID<DIMENSION> gid;
	gid.node_idx = nid;
	gid.direction = dir;
	typename GhostMap::iterator it = ghostmap.find(gid);
	return it->second;
}

template<class DIMENSION>
void BCManager<DIMENSION>::set_bc_on_ghost_nodes(
		typename BCManager<DIMENSION>::GhostMap::iterator& it_gn,
		typename BCManager<DIMENSION>::BCSet::iterator& it_bc) {
	it_bc->pfun((*it_bc), it_gn->second);
}
template<class DIMENSION>
void BCManager<DIMENSION>::set_bc() {
	for (typename GhostMap::iterator it_gn = ghostmap.begin();
			it_gn != ghostmap.end(); ++it_gn) {
		for (typename BCSet::iterator it_bc = bcset.begin();
				it_bc != bcset.end(); ++it_bc) {
			if (it_gn->first.tree_idx == it_bc->tree_idx
					&& it_gn->first.direction == it_bc->direction) {
				set_bc_on_ghost_nodes(it_gn, it_bc);
			}
		}
	}
}

} //end namespace

#endif /* CALCULATION_BOUNDARY_H_ */
