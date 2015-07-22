/*
 * Advection.cpp
 *
 *  Created on: Mar 11, 2015
 *      Author: zhou
 */

#include "Advection.h"
#include "Scalar.h"
#include "../Utility/ArrayList.h"
#include "../IO/IO.h"
#include "../IO/IO_vtk.h"
#include "../Algebra/Arithmetic.h"
#include "../Algebra/Interpolation.h"

namespace Larus {

inline Float _gernal_k_scheme(Float c, Float dc, Float gp, Float gm, Float k) {
	return c + dc / 2.0 * ((1 + k) / 2.0 * gp + (1 - k) / 2.0 * gm);
}

Float k_scheme_CDS(Float c, Float dc, Float gp, Float gm) { //central difference
	return _gernal_k_scheme(c, dc, gp, gm, 1.0);
}
Float k_scheme_QUICK(Float c, Float dc, Float gp, Float gm) { //Quadratic-upwind
	return _gernal_k_scheme(c, dc, gp, gm, 0.5);
}
Float k_scheme_CUI(Float c, Float dc, Float gp, Float gm) {    //cubic-upwind
	return _gernal_k_scheme(c, dc, gp, gm, 1.0 / 3.0);
}
Float k_scheme_Fromm(Float c, Float dc, Float gp, Float gm) {  //Fromm's scheme
	return _gernal_k_scheme(c, dc, gp, gm, 0.0);
}
Float k_scheme_SOU(Float c, Float dc, Float gp, Float gm) { //second-order upwind
	return _gernal_k_scheme(c, dc, gp, gm, -1.0);
}

//===============================================

void _subvisit_traversal_to_ListFace(pQTNode pn, pQTNode pnei,
		ListT<QTNodeFace>& listf, SPDirection d, SPNodeFaceType ft) {

	if (ft == SPFT_FineCoarse || ft == SPFT_Boundary) {
		QTNodeFace nf(pn, pnei, d, ft);
		listf.push_back(nf);
	} else { //ft == SPFT_Equal
		if (d == SPD_IM || d == SPD_IP) {
			Float veof = interpolate_1order_on_face(pn, d, U_IDX);
			if (d == SPD_IP && veof >= 0) {
				QTNodeFace nf(pn, pnei, d, ft);
				listf.push_back(nf);
			}
			if (d == SPD_IM && veof <= 0) {
				QTNodeFace nf(pn, pnei, d, ft);
				listf.push_back(nf);
			}
		}
		if (d == SPD_JM || d == SPD_JP) {
			Float veof = interpolate_1order_on_face(pn, d, V_IDX);
			if (d == SPD_JP && veof >= 0) {
				QTNodeFace nf(pn, pnei, d, ft);
				listf.push_back(nf);
			}
			if (d == SPD_JM && veof <= 0) {
				QTNodeFace nf(pn, pnei, d, ft);
				listf.push_back(nf);
			}
		}
	}
}

void _visit_traversal_to_ListFace(pQTNode pn, utPointer utp) {
	if (condition_is_leaf(pn)) {
		for (int i = 4; i <= 7; ++i) {
			SPDirection dir = toDirection(i);
			pQTNode pnei = pn->getNeighborFast(dir);
			SPNodeFaceType ft = getFaceType(pn, pnei);
			if(ft == SPFT_Error || ft == SPFT_CoarseFine)
				continue;
			ListT<QTNodeFace>& listf = *CAST(ListT<QTNodeFace>*, utp);
			_subvisit_traversal_to_ListFace(pn, pnei, listf, dir, ft);
		}
	}
}

void traversal_to_ListFace(pQuadTree ptree, ListT<QTNodeFace>& listf) {
	_IF_TRUE_RETRUN(ptree==NULL_PTR);
	ptree->Traversal(_visit_traversal_to_ListFace, &listf);
}

void draw_gnuplot_ListFace(std::string filename,  //filename
		ListT<QTNodeFace>& listf, //list
		int mode              //mode
		) {
	FILE *data = open_file(filename, mode);
	for (ListT<QTNodeFace>::iterator iter = listf.begin(); iter != listf.end();
			++iter) {
		if (iter->pnode != NULL_PTR && iter->face_type != SPFT_Error) {
			if (iter->direction == SPD_IM) {
				fprintf(data, "%lf %lf \n",
						iter->pnode->cell->getPoint(eCPL_M, eCPL_P).x,
						iter->pnode->cell->getPoint(eCPL_M, eCPL_P).y);
				fprintf(data, "%lf %lf \n",
						iter->pnode->cell->getPoint(eCPL_M, eCPL_M).x,
						iter->pnode->cell->getPoint(eCPL_M, eCPL_M).y);
			}
			if (iter->direction == SPD_IP) {
				fprintf(data, "%lf %lf \n",
						iter->pnode->cell->getPoint(eCPL_P, eCPL_P).x,
						iter->pnode->cell->getPoint(eCPL_P, eCPL_P).y);
				fprintf(data, "%lf %lf \n",
						iter->pnode->cell->getPoint(eCPL_P, eCPL_M).x,
						iter->pnode->cell->getPoint(eCPL_P, eCPL_M).y);
			}
			if (iter->direction == SPD_JM) {
				fprintf(data, "%lf %lf \n",
						iter->pnode->cell->getPoint(eCPL_M, eCPL_M).x,
						iter->pnode->cell->getPoint(eCPL_M, eCPL_M).y);
				fprintf(data, "%lf %lf \n",
						iter->pnode->cell->getPoint(eCPL_P, eCPL_M).x,
						iter->pnode->cell->getPoint(eCPL_P, eCPL_M).y);
			}
			if (iter->direction == SPD_JP) {
				fprintf(data, "%lf %lf \n",
						iter->pnode->cell->getPoint(eCPL_M, eCPL_P).x,
						iter->pnode->cell->getPoint(eCPL_M, eCPL_P).y);
				fprintf(data, "%lf %lf \n",
						iter->pnode->cell->getPoint(eCPL_P, eCPL_P).x,
						iter->pnode->cell->getPoint(eCPL_P, eCPL_P).y);
			}
			fprintf(data, "\n");
		}
	}

	fclose(data);
}

}

