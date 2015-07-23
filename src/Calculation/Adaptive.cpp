/*
 * Adaptive.cpp
 *
 *  Created on: Jul 6, 2015
 *      Author: zhou
 */

#include "Adaptive.h"

namespace Larus {

inline bool level_check(int minlevel, int maxlevel) {
	// the max level is great than the min level
	// and the min level is great than 1
	if (minlevel <= maxlevel && minlevel >= 1) {
		return true;
	} else {
		return false;
	}
}
template<class FOREST>
void creat_forest_to_level(FOREST& f, int level) {
	typedef typename FOREST::Tree Tree;
	for (int i = 0; i < f.size(); i++) {
		Tree* tree = f.getpTree_1d(i);
		if (tree != NULL_PTR) {
			tree->CreatFullTree_to_Level(level);
		}
	}
}
template<class NODE>
void visit_refine_point_to_level(NODE* pn, utPointer p) {
	typedef SPTree<NODE, NODE::DIM>* ptree_t;
	arrayListT<utPointer>& arrp = (*CAST(arrayListT<utPointer>*, p));
	ptree_t pt = CAST(ptree_t, arrp[0]);
	int& level = (*CAST(int*, arrp[1]));
	typedef typename SPTree<NODE, NODE::DIM>::Node::Cell_type::Point Point;
	Point& poi = (*CAST(Point*, arrp[2]));
	if (condition_is_node(pn) && pn->getLevel() < level
			&& pn->cell->isInOnCell(poi)) {
		pt->Creatchildren_full(pn);
	}
}

void initial_adaptation(Forest2D& forest,      // Forest 2D
		Point2D& point,        // A Point2D
		int minlevel,          // minlevel
		int maxlevel) {         // maxlevel
	assert(level_check(minlevel, maxlevel));
	// we are assume that the forest did not been initialized
	// 1 create tree to minlevel
	creat_forest_to_level(forest, minlevel);
	// 2 refine the location where the point is located.
	for (int i = 0; i < forest.size(); i++) {
		pQuadTree tree = forest.getpTree_1d(i);
		if (tree != NULL_PTR) {
			if (tree->isInOnTree(point)) {
				arrayListT<utPointer> arr(3);
				arr[0] = tree;
				arr[1] = &maxlevel;
				arr[2] = &point;
				tree->Traversal(visit_refine_point_to_level, &arr);
			}
		}
	}

}

void initial_adaptation(Forest2D& forest,      // Forest 2D
		ListT<Point2D>& lp,        // A Point2D
		int minlevel,          // minlevel
		int maxlevel) {         // maxlevel
	assert(level_check(minlevel, maxlevel));
	// we are assume that the forest did not been initialized
	// 1 create tree to minlevel
	creat_forest_to_level(forest, minlevel);
	// 2 refine the location where the point is located.
	for (int i = 0; i < forest.size(); i++) {
		pQuadTree tree = forest.getpTree_1d(i);
		if (tree != NULL_PTR) {
			for (auto iter = lp.begin(); iter != lp.end(); ++iter) {
				if (tree->isInOnTree((*iter))) {
					arrayListT<utPointer> arr(3);
					Point2D* pp = &(*iter);
					arr[0] = tree;
					arr[1] = &maxlevel;
					arr[2] = pp;
					tree->Traversal(visit_refine_point_to_level, &arr);
				}
			}
		}
	}
}

template<class CELL>
bool _is_refine_eq(CELL& cell, pFun pf, Float value) {
		int flag = sign(
				value
						- pf(cell.getX(CELL_POINT_ORDER[0][0]),
								cell.getY(CELL_POINT_ORDER[0][1]),
								cell.getZ(CELL_POINT_ORDER[0][2])));
		for (int i = 1; i < CELL::NUM_VERTEXES; i++) {
			int b = sign(
					value
							- pf(cell.getX(CELL_POINT_ORDER[i][0]),
									cell.getY(CELL_POINT_ORDER[i][1]),
									cell.getZ(CELL_POINT_ORDER[i][2])));
			if (flag != b) {
				return true;
			}
		}
		return false;

}

template<class CELL>
bool _is_refine_le(CELL& cell, pFun pf, Float value) { //less equal
	for (int i = 1; i < CELL::NUM_VERTEXES; i++) {
		int flag = sign(
				value
						- pf(cell.getX(CELL_POINT_ORDER[i][0]),
								cell.getY(CELL_POINT_ORDER[i][1]),
								cell.getZ(CELL_POINT_ORDER[i][2])));
		if (flag == 1) {
			return true;
		}
	}
	return false;
}

template<class NODE>
void visit_refine_function_eq(NODE* pn, utPointer p) {
	typedef SPTree<NODE, NODE::DIM>* ptree_t;
	arrayListT<utPointer>& arrp = (*CAST(arrayListT<utPointer>*, p));
	ptree_t pt = CAST(ptree_t, arrp[0]);
	int& level = (*CAST(int*, arrp[1]));
	pFun& pf = (*CAST(pFun*, arrp[2]));
	Float& value = (*CAST(Float*, arrp[3]));
	if (condition_is_node(pn) && pn->getLevel() < level
			&& _is_refine_eq(*(pn->cell), pf, value)) {
		pt->Creatchildren_full(pn);
	}
}

template<class NODE>
void visit_refine_function_le(NODE* pn, utPointer p) {
	typedef SPTree<NODE, NODE::DIM>* ptree_t;
	arrayListT<utPointer>& arrp = (*CAST(arrayListT<utPointer>*, p));
	ptree_t pt = CAST(ptree_t, arrp[0]);
	int& level = (*CAST(int*, arrp[1]));
	pFun& pf = (*CAST(pFun*, arrp[2]));
	Float& value = (*CAST(Float*, arrp[3]));
	if (condition_is_node(pn) && pn->getLevel() < level
			&& _is_refine_le(*(pn->cell), pf, value)) {
		pt->Creatchildren_full(pn);
	}
}

void initial_adaptation_eq( //
		Forest2D& forest,      // Forest 2D
		pFun pf,               // A Point2D
		Float value,           // Threshold value
		int minlevel,          // minlevel
		int maxlevel) {         // maxlevel
	assert(level_check(minlevel, maxlevel));
	// we are assume that the forest did not been initialized
	// 1 create tree to minlevel
	creat_forest_to_level(forest, minlevel);
	// 2 refine the location where the function is satisfied
	for (int i = 0; i < forest.size(); i++) {
		pQuadTree tree = forest.getpTree_1d(i);
		if (tree != NULL_PTR) {
			//if (_is_refine_eq(*(tree->getpRootCell()), pf, value)) {
				arrayListT<utPointer> arr(4);
				arr[0] = tree;
				arr[1] = &maxlevel;
				arr[2] = &pf;
				arr[3] = &value;
				tree->Traversal(visit_refine_function_eq, &arr);
			//}
		}

	}
}

void initial_adaptation_le( //
		Forest2D& forest,      // Forest 2D
		pFun pf,               // A Point2D
		Float value,           // Threshold value
		int minlevel,          // minlevel
		int maxlevel) {         // maxlevel
	assert(level_check(minlevel, maxlevel));
	// we are assume that the forest did not been initialized
	// 1 create tree to minlevel
	creat_forest_to_level(forest, minlevel);
	// 2 refine the location where the function is satisfied
	for (int i = 0; i < forest.size(); i++) {
		pQuadTree tree = forest.getpTree_1d(i);
		if (tree != NULL_PTR) {
			//if (_is_refine_eq(*(tree->getpRootCell()), pf, value)) {
				arrayListT<utPointer> arr(4);
				arr[0] = tree;
				arr[1] = &maxlevel;
				arr[2] = &pf;
				arr[3] = &value;
				tree->Traversal(visit_refine_function_le, &arr);
			//}
		}

	}
}

}
