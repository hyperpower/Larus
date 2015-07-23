/************************
 //  \file   Levelset.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   10 févr. 2015 
 ***********************/

#include "../TypeDef.h"
#include "Levelset.h"
#include "../Geometry/Relation.h"

namespace Larus {
//work with tree=================================
//===============================================
//set value======================================

void visit_set_sphere_ls(OCNode* node, utPointer pt) {
	if (condition_is_leaf(node)) {
		arrayListT<utPointer> arr = (*CAST(arrayListT<utPointer>*, pt));
		Point3D& p = (*CAST(Point3D*, arr[0]));
		Float& r = (*CAST(Float*, arr[1]));
		LarusDef::size_type& idx = (*CAST(LarusDef::size_type*, arr[2]));
		node->data->aCenterData[idx] = calDistance(node->cell->getCenterPoint(),
				p) - r;
	}
}
void set_sphere_ls(pOCTree oct, Point3D p, Float r, LarusDef::size_type idx) {
	arrayListT<utPointer> arr(3);
	arr[0] = &p;
	arr[1] = &r;
	arr[2] = &idx;
	oct->Traversal(visit_set_sphere_ls, &arr);
}

bool is_sphere_boundary_in_cell(pOCNode pnode, const Point3D& p, Float r) {
	Float dis = ABS(calDistance(pnode->cell->getCenterPoint(), p) - r);
	if (dis < pnode->cell->gethDx() || dis < pnode->cell->gethDy()
			|| dis < pnode->cell->gethDz()) {
		return true;
	} else if (dis > pnode->cell->gethDiagnalLength()) {
		return false;
	} else {
		bool bm = is_in(pnode->cell->getVertex(0), p, r);
		for (LarusDef::size_type i = 1; i < pnode->cell->NUM_VERTEXES; i++) {
			if (bm != is_in(pnode->cell->getVertex(i), p, r)) {
				return true;
			}
		}
		return false;
	}
}

int condition_re_sphere_ls(arrayList& avt, pOCNode pnode, utPointer pt) {
	arrayListT<utPointer> arr = (*CAST(arrayListT<utPointer>*, pt));
	Point3D& p = (*CAST(Point3D*, arr[0]));
	Float& r = (*CAST(Float*, arr[1]));
	//LarusDef::size_type& idx = (*CAST(LarusDef::size_type*, arr[2]));
	if (is_sphere_boundary_in_cell(pnode, p, r)) {
		avt.assign(1);
		return 1;
	} else {
		return -1;
	}
}

void set_sphere_ls_initial_refine(pOCTree oct, Point3D p, Float r,
		LarusDef::size_type idx) {
	arrayListT<utPointer> arr(3);
	arr[0] = &p;
	arr[1] = &r;
	arr[2] = &idx;
	oct->Refine(condition_re_sphere_ls, &arr, visit_empty, &arr);
}

}

