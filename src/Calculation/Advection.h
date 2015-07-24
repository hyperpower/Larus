/*
 * Advection.h
 *
 *  Created on: Mar 11, 2015
 *      Author: zhou
 */

#ifndef _ADVECTION_H_
#define _ADVECTION_H_

#include "../TypeDef.h"
#include "CalDef.h"
#include "../Grid/SPTree.h"
#include "../Grid/SPTreeNode.h"
#include "../Grid/Forest.h"
#include "../Algebra/Arithmetic.h"
#include "Exp.h"
#include "Boundary.h"

namespace Larus {

//This file use to solve advection equation
//
//   d(phi)      d(phi)      d(phi)       2D --version
//   ------ + u ------- + v -------  = 0
//     dt         dx           dy
//
//   d(phi)     d(phi)       d(phi)       d(phi)       3D --version
//   ------ + u ------- + v -------  + w -------  = 0
//     dt         dx           dy          dz

template<class DIMENSION>
class Advection_Eq: public ObjectBase {
public:
	typedef typename DIMENSION::value_type vt;
	typedef typename DIMENSION::size_type st;
	typedef Forest<DIMENSION> Forest_;
	typedef typename Forest_::pFace pFace;
	typedef typename Forest_::Face Face;
	typedef typename DIMENSION::pTree pTree;
	typedef typename DIMENSION::Tree Tree;
	typedef typename DIMENSION::pNode pNode;
	typedef typename DIMENSION::Node Node;
	typedef typename DIMENSION::Point Point;
	typedef ExpTNode<vt, Node> ExpTerm;
	typedef Exp<ExpTerm> Expression;

	typedef vt (*pfun_f)(vt, vt, vt);

	Forest_* pforest;
	BCManager<DIMENSION>* pBCM;

	Float t;

	st phi_idx;
	st dt_idx;
	st u_idx;
	st v_idx;
	st w_idx;

	Advection_Eq(Forest_*, BCManager<DIMENSION>*, st, st, st, st, st = 0);

	// set --------------------------------------
	void _set_val(st, pfun_f);
	void set_phi(pfun_f pfun) {
		_set_val(phi_idx, pfun);
	}
	void set_dt(pfun_f pfun) {
		_set_val(dt_idx, pfun);
	}
	void set_u(pfun_f pfun) {
		_set_val(u_idx, pfun);
	}
	void set_v(pfun_f pfun) {
		_set_val(v_idx, pfun);
	}
	void set_w(pfun_f pfun) {
		ASSERT(DIMENSION::DIM == 3);
		_set_val(w_idx, pfun);
	}

	int _face_scheme_equal_adv(pFace, Expression&);
	int _face_scheme_adv(pFace, Expression&);
	int _face_exp_adv(pNode);

};

template<class DIMENSION>
Advection_Eq<DIMENSION>::Advection_Eq(Forest_* pf, BCManager<DIMENSION>* pBCM,
		st phii, st ti, st ui, st vi, st wi) :
		pforest(pf), pBCM(pBCM), phi_idx(phii), dt_idx(ti), u_idx(ui), v_idx(
				vi), w_idx(wi) {
	t = 0; //initial time to 0.0
}

template<class DIMENSION>
void Advection_Eq<DIMENSION>::_set_val(st idx, pfun_f pfun) {
	for (typename Forest_::iterator it = pforest->begin(); it != pforest->end();
			++it) {
		if (DIMENSION::DIM == 2) {
			it->data->aCenterData[idx] = pfun(it->cell->getCPX(),
					it->cell->getCPY(), 0);
		} else {
			it->data->aCenterData[idx] = pfun(it->cell->getCPX(),
					it->cell->getCPY(), it->cell->getCPZ());
		}
	}
}
template<class DIMENSION>
int Advection_Eq<DIMENSION>::_face_exp_adv(pNode pn) {
	typedef Pair<Face*, Expression*> Pair;
	typedef ListT<Pair> List;
	if (DIMENSION::DIM == 2) {
		if (pn->data->utp_data == NULL_PTR) {
			pn->data->utp_data = new List();   //new !!!!!
		}
		for (int i = 4; i <= 7; i++) {
			SPDirection d = toDirection(i);
			pNode pnei = pn->getNeighborFast(d);
			// The face need to be calculate
			// 1 the boundary face
			// 2 the equal face on P direction;
			// 3 the fine to coarse face

			// case 1 2 3
			Face* f = new Face(pn, pnei, d, getFaceType(pn, pnei));
			Expression* pexp = new Expression();
			if ((f->face_type == SPFT_Boundary)
					|| (f->face_type == SPFT_Equal
							&& ((d == SPD_IP) || (d == SPD_JP))) //2
					|| (f->face_type == SPFT_FineCoarse))  //1
					{
				// work on pn
				_face_scheme_adv(f, *pexp);
				List& lpexp = (*CAST(List*, pn->data->utp_data));
				Pair pair(f, pexp);
				lpexp.push_back(pair);

			}
			// case 2 3
			if ((f->face_type == SPFT_Equal && ((d == SPD_IP) || (d == SPD_JP))) //2
			|| (f->face_type == SPFT_FineCoarse)) //3
					{
				// work on pnei
				if (pnei->data->utp_data == NULL_PTR) {
					pnei->data->utp_data = new List();    //new !!!!!
				}
				SPNodeFaceType ft =
						(f->face_type == SPFT_FineCoarse) ?
								SPFT_CoarseFine : SPFT_Equal;
				Face* fn = new Face(pnei, pn, oppositeDirection(d), ft);
				Expression* pexpn = new Expression(*pexp);
				//pexpn->times(-1.0);
				List& lpexp = (*CAST(List*, pnei->data->utp_data));
				Pair pairn(fn, pexpn);
				lpexp.push_back(pairn);
			}
		}
	}
	return 1;
}

template<class DIMENSION>
int Advection_Eq<DIMENSION>::_face_scheme_adv(pFace pface, Expression& exp) {
	//face type
	switch (pface->face_type) {
	case SPFT_Error:
		break;
	case SPFT_Boundary:
		std::cout << "Boundary function unfinish\n";
		break;
	case SPFT_Equal: {
		this->_face_scheme_equal_adv(pface, exp);
		break;
	}
	case SPFT_FineCoarse:
		std::cout << "Fine Corase function unfinish\n";
		break;
	case SPFT_CoarseFine:
		std::cout << "Corase Fine function unfinish\n";
		break;
	default:
		return -1;
	}
	return pface->face_type;
}
// This function calculate the phi value on face
template<class DIMENSION>
int Advection_Eq<DIMENSION>::_face_scheme_equal_adv(pFace pface,
		Expression& exp) {
	// face direction in 4,5,6,7
	//                  4         5         6         7
	CSAxis arr_dd[] = { CSAxis_X, CSAxis_Y, CSAxis_X, CSAxis_Y, CSAxis_Z,
			CSAxis_Z };
	int arr_signd[] = { -1, 1, 1, -1, 1, -1 };
	int arr_veo_idx[] = { u_idx, v_idx, u_idx, v_idx, w_idx, w_idx };

	CSAxis dd = arr_dd[int(pface->direction) - 4];
	st veo_idx = arr_veo_idx[int(pface->direction) - 4];
	Float veo_f = interpolate_1order_on_face((*pface), veo_idx);
	//get U C D ------------------------------
	pNode pU=NULL_PTR;
	pNode pC=NULL_PTR;
	pNode pD=NULL_PTR;
	if(isPlus(pface->direction)){
		if(veo_f>0){
			pC = pface->node;
			pD = pface->neighbor;
		}else{

		}
	}else{
		if(veo_f>0){

		}else{

		}
	}
}

//k-scheme ======================================
//===============================================
Float k_scheme_CDS(Float c, Float dc, Float gp, Float gm);  //central difference
Float k_scheme_QUICK(Float c, Float dc, Float gp, Float gm);  //Quadratic-upwind
Float k_scheme_CUI(Float c, Float dc, Float gp, Float gm);    //cubic-upwind
Float k_scheme_Fromm(Float c, Float dc, Float gp, Float gm);  //Fromm's scheme
Float k_scheme_SOU(Float c, Float dc, Float gp, Float gm); //second-order upwind

//Traversal ListFace=============================
void traversal_to_ListFace(pQuadTree ptree, ListT<QTNodeFace>& listf);
//draw ==========================================
void draw_gnuplot_ListFace(std::string filename,  //filename
		ListT<QTNodeFace>& listf, //list
		int mode              //mode
		);

} //end of namespace

#endif /* CALULATION_ADVECTION_H_ */
