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
#include "../Algebra/MatrixSparCompRow.h"
#include "../Algebra/Solver_matrix.h"
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

typedef Float (*pfun_limiter)(Float);

// First order up wind
inline Float FOU(Float r) {
	return 0;
}
inline Float limiter_minmod(Float r) {
	return MAX(0.0, MIN(1.0, r));  // unfinish
}
inline Float limiter_superbee(Float r) {
	return MAX(0.0, MIN(2.0 * r, 1.0), MIN(r, 2.0));
}
inline Float limiter_vanLeer(Float r) {
	return (r + ABS(r)) / (1 + r);
}

const pfun_limiter LIMITER_LIST[] = {  //
		FOU,  //0
				limiter_minmod,  //1
				limiter_superbee, //2
				limiter_vanLeer };

inline Float limiter(Float r, int i) {
	return LIMITER_LIST[i](r);
}

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
	typedef BCManager<DIMENSION> BCM;
	typedef ExpTNode<vt, Node> ExpTerm;
	typedef Exp<ExpTerm> Expression;

	typedef vt (*pfun_f)(vt, vt, vt);

	Forest_* pforest;
	BCM* pBCM;

	Float t;
	st scheme_idx;
	st phi_idx;
	st phin_idx;
	st phinn_idx;
	st dt_idx;
	st u_idx;
	st v_idx;
	st w_idx;

	Advection_Eq(Forest_* pf, BCManager<DIMENSION>* pBCM, st scheme, //
			st phii, st phini, st phinni, //
			st dti, st ui, st vi, st wi = 0) :
			pforest(pf), //
			pBCM(pBCM),  //
			scheme_idx(scheme), //
			phi_idx(phii), //
			phin_idx(phini), //
			phinn_idx(phinni), //
			dt_idx(dti), //
			u_idx(ui), //
			v_idx(vi), //
			w_idx(wi) {
		t = 0; //initial time to 0.0
	}

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

	int _substitute_boudary_val(Expression& exp, st idx);

	// first order up wind ------------------
	int _find_C(pFace, pNode&, Float, BCM*);
	int _find_C(pFace, pNode&, Float);
	int _face_scheme_boundary_fou(pFace, Expression&);
	int _face_scheme_equal_fou(pFace, Expression&);
	int _face_scheme_fine_coarse_fou(pFace, Expression&);
	int _face_scheme_fou(pFace, Expression&);
	int _face_exp_fou(pNode);

	Float _face_scheme_boundary_tvd(pFace, Expression&, st);
	Float _face_scheme_equal_tvd(pFace, Expression&, st);
	Float _face_scheme_fine_coarse_tvd(pFace, Expression&, st);
	Float _face_scheme_tvd(pFace, Expression&, st);
	int _face_exp_tvd(pNode, st);

	int _node_exp_adv_t(pNode, Float, Expression&);
	int _node_exp_adv(pNode, Expression&);

	int _clear_utp_data(pNode);
	int advance(int step = 1);

	int _bulid_matrix_fou(MatrixSCR<Float>& mat, arrayListV<Float>& b);
	int _bulid_matrix_tvd(MatrixSCR<Float>& mat, arrayListV<Float>& b);
	int slove(Float tol);
	int slove_fou(Float tol);
};

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
int Advection_Eq<DIMENSION>::_face_exp_fou(pNode pn) {
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
			Face* pf = new Face(pn, pnei, d, getFaceType(pn, pnei));
			Expression* pexp = new Expression();
			if ((pf->face_type == SPFT_Boundary)            //1
			|| (pf->face_type == SPFT_Equal && ((d == SPD_IP) || (d == SPD_JP))) //2
					|| (pf->face_type == SPFT_FineCoarse))  //3
					{
				// work on pn
				_face_scheme_fou(pf, *pexp);
				List& lpexp = (*CAST(List*, pn->data->utp_data));
				Pair pair(pf, pexp);
				lpexp.push_back(pair);
			}
			// case 2 3
			if ((pf->face_type == SPFT_Equal && ((d == SPD_IP) || (d == SPD_JP))) //2
			|| (pf->face_type == SPFT_FineCoarse)) //3
					{
				// work on pnei
				if (pnei->data->utp_data == NULL_PTR) {
					pnei->data->utp_data = new List();    //new !!!!!
				}
				SPNodeFaceType ft =
						(pf->face_type == SPFT_FineCoarse) ?
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
int Advection_Eq<DIMENSION>::_face_exp_tvd(pNode pn, st idx) {
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
				_face_scheme_tvd(f, *pexp, idx);
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
Float Advection_Eq<DIMENSION>::_face_scheme_tvd(pFace pface, Expression& exp,
		st idx) {
	//face type
	switch (pface->face_type) {
	case SPFT_Error:
		break;
	case SPFT_Boundary:
		this->_face_scheme_boundary_tvd(pface, exp, idx);
		break;
	case SPFT_Equal: {
		this->_face_scheme_equal_tvd(pface, exp, idx);
		break;
	}
	case SPFT_FineCoarse:
		this->_face_scheme_fine_coarse_tvd(pface, exp, idx);
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
Float Advection_Eq<DIMENSION>::_face_scheme_equal_tvd(pFace pface,
		Expression& exp, st idx) {
	int arr_veo_idx[] = { u_idx, v_idx, u_idx, v_idx, w_idx, w_idx };

	st veo_idx = arr_veo_idx[int(pface->direction) - 4];
	Float veo_f = interpolate_1order_on_face((*pface), veo_idx);
	//get U C D ------------------------------
	pNode pU = NULL_PTR;
	pNode pC = NULL_PTR;
	pNode pD = NULL_PTR;
	if (isPlus(pface->direction)) {
		if (veo_f > 0) {
			pC = pface->pnode;
			pD = pface->pneighbor;
			find_neighbor_2(pC, oppositeDirection(pface->direction), pU,
					(*pBCM));
		} else {
			pC = pface->pneighbor;
			pD = pface->pnode;
			find_neighbor_2(pC, pface->direction, pU, (*pBCM));
		}
	} else {
		if (veo_f > 0) {
			pC = pface->pneighbor;
			pD = pface->pnode;
			find_neighbor_2(pC, pface->direction, pU, (*pBCM));
		} else {
			pC = pface->pnode;
			pD = pface->pneighbor;
			find_neighbor_2(pC, oppositeDirection(pface->direction), pU,
					(*pBCM));
		}
	}
	//
	// exp should be empty
	exp.Insert(ExpTerm(getIDX(pC), pC, 1.0));
	// Analyze pU
	//Expression expU;
	//assert(pU != NULL_PTR);
	//if (condition_is_leaf(pU)) {
	//	expU.Insert(ExpTerm(getIDX(pU), pU, 1.0));  //
	//} else { //pU has children -------
	//	ListT<pQTNode> lch;
	//	getListpNode_children(pU, lch);
	//	ListT<pQTNode>::size_type count = lch.size();
	//	for(ListT<pQTNode>::iterator iter=lch.begin(); iter!=lch.end(); ++iter){
	//		pQTNode pn=(*iter);
	//		expU.Insert(ExpTerm(getIDX(pn), pn, 1.0/double(count)));
	//	}
	//}
	//--------------------------
	// cal limiter
	Float vU = getAverageVal(pU, idx);
	Float vC = getcVal(pC, idx);
	Float vD = getcVal(pD, idx);
	// cal \Psi(r)
	Float r = (vC - vU) / (vD - vC + SMALL);
	Float psi = limiter(r, scheme_idx);
	Float cor = 0.5 * psi * (vD - vC);
	exp.Insert(ExpTerm(ExpTerm::IDX_CONST, NULL_PTR, cor));

	return cor;
}

template<class DIMENSION>
Float Advection_Eq<DIMENSION>::_face_scheme_boundary_tvd(pFace pface,
		Expression& exp, st idx) {
	// face direction in 4,5,6,7
	//                  4         5         6         7
	//CSAxis arr_dd[] = { CSAxis_X, CSAxis_Y, CSAxis_X, CSAxis_Y, CSAxis_Z,
	//		CSAxis_Z };
	//int arr_signd[] = { -1, 1, 1, -1, 1, -1 };
	int arr_veo_idx[] = { u_idx, v_idx, u_idx, v_idx, w_idx, w_idx };

	//CSAxis dd = arr_dd[int(pface->direction) - 4];
	st veo_idx = arr_veo_idx[int(pface->direction) - 4];
	Float veo_f = interpolate_1order_on_face((*pface), veo_idx);
	//get U C D ------------------------------
	pNode pU = NULL_PTR;
	pNode pC = NULL_PTR;
	pNode pD = NULL_PTR;
	if (isPlus(pface->direction)) {
		if (veo_f > 0) {
			pC = pface->pnode;  //o
			pD = pBCM->find_ghost(getIDX(pface->pnode), pface->direction);
			find_neighbor_2(pC, oppositeDirection(pface->direction), pU,
					(*pBCM));
		} else {
			pC = pBCM->find_ghost(getIDX(pface->pnode), pface->direction);
			pD = pface->pnode;
			pU = pC;
		}
	} else {
		if (veo_f > 0) {
			pC = pBCM->find_ghost(getIDX(pface->pnode), pface->direction);
			pD = pface->pnode;
			pU = pC;
		} else {
			pC = pface->pnode;
			pD = pBCM->find_ghost(getIDX(pface->pnode), pface->direction);
			find_neighbor_2(pC, oppositeDirection(pface->direction), pU,
					(*pBCM));
		}
	}
	// exp should be empty
	exp.Insert(ExpTerm(getIDX(pC), pC, 1.0));
	// Analyze pU
	ASSERT(pU!= NULL_PTR);
	// cal limiter
	Float vU = getAverageVal(pU, idx);
	Float vC = getcVal(pC, phi_idx);
	Float vD = getcVal(pD, phi_idx);
	// cal \Psi(r)
	Float r = (vC - vU) / (vD - vC + SMALL);
	Float psi = limiter(r, scheme_idx);
	Float cor = 0.5 * psi * (vD - vC);
	exp.Insert(ExpTerm(ExpTerm::IDX_CONST, NULL_PTR, cor));

	return cor;

}

template<class DIMENSION>
Float Advection_Eq<DIMENSION>::_face_scheme_fine_coarse_tvd(pFace pface,
		Expression& exp, st idx) {
	int arr_veo_idx[] = { u_idx, v_idx, u_idx, v_idx, w_idx, w_idx };

	st veo_idx = arr_veo_idx[int(pface->direction) - 4];
	Float veo_f = interpolate_1order_on_face((*pface), veo_idx);
	//get U C D ------------------------------
	pNode pU = NULL_PTR;
	pNode pC = NULL_PTR;
	pNode pD = NULL_PTR;
	if (isPlus(pface->direction)) {
		if (veo_f > 0) {
			pC = pface->pnode;
			pD = pface->pneighbor;
			find_neighbor_2(pC, oppositeDirection(pface->direction), pU,
					(*pBCM));
		} else {
			pC = pface->pneighbor;
			pD = pface->pnode;
			find_neighbor_2(pC, pface->direction, pU, (*pBCM));
		}
	} else {
		if (veo_f > 0) {
			pC = pface->pneighbor;
			pD = pface->pnode;
			find_neighbor_2(pC, pface->direction, pU, (*pBCM));
		} else {
			pC = pface->pnode;
			pD = pface->pneighbor;
			find_neighbor_2(pC, oppositeDirection(pface->direction), pU,
					(*pBCM));
		}
	}
	//
	// exp should be empty
	exp.Insert(ExpTerm(getIDX(pC), pC, 1.0));

	// cal limiter
	Float vU = getAverageVal(pU, idx);
	Float vC = getcVal(pC, idx);
	Float vD = getcVal(pD, idx);
	// cal \Psi(r)
	Float r = (vC - vU) / (vD - vC + SMALL);
	Float psi = limiter(r, scheme_idx);
	Float cor = 0.5 * psi * (vD - vC);
	exp.Insert(ExpTerm(ExpTerm::IDX_CONST, NULL_PTR, cor));

	return cor;
}


template<class DIMENSION>
int Advection_Eq<DIMENSION>::_node_exp_adv_t(pNode pn, Float dt,
		Expression& exp) {
	//typedef typename Forest_::Node Node;
	typedef typename Forest_::Face Face;
	typedef Pair<Face*, Expression*> Pair;
	typedef ListT<Pair> List;

	if (DIMENSION::DIM == 2) {
		assert(pn->data->utp_data !=NULL_PTR);
		List& lpexp = (*CAST(List*, pn->data->utp_data));
		Expression FF[4];
		int countFF[] = { 0, 0, 0, 0 };
		for (typename List::iterator iter = lpexp.begin(); iter != lpexp.end();
				++iter) {
			assert(iter->first->pnode == pn); //
			Face* iterf = iter->first;
			Expression* itere = iter->second;
			int idd = int(iterf->direction) - 4;
			FF[idd].plus(*itere);
			countFF[idd]++;
		}
		for (int i = 0; i < 4; i++) {
			if (countFF[i] > 1) {
				FF[i].times(1.0 / Float(countFF[i])); // get average face gradient;
			}
		}
		// x direction
		Expression exp_x = FF[2] - FF[0];
		//std::cout<<" exp _ x" << exp_x.size()<<"\n";
		Float l = pn->cell->getD(CSAxis_X);
		Float u = getcVal(pn, u_idx);
		Float CFL_x = u * dt / l;
		ASSERT_MSG(CFL_x < 0.5, " CFL x > 0.5");
		exp_x.times(-CFL_x);  //negative
		// y direction
		Expression exp_y = FF[1] - FF[3];
		l = pn->cell->getD(CSAxis_Y);
		u = getcVal(pn, v_idx);
		Float CFL_y = u * dt / l;
		ASSERT_MSG(CFL_y < 0.5, " CFL y > 0.5");
		exp_y.times(-CFL_y); //negative
		//
		exp.plus(exp_x);
		exp.plus(exp_y);

		_clear_utp_data(pn);
	}
	return 1;
}

template<class DIMENSION>
int Advection_Eq<DIMENSION>::_node_exp_adv(pNode pn, Expression& exp) {
	//typedef typename Forest_::Node Node;
	typedef typename Forest_::Face Face;
	typedef Pair<Face*, Expression*> Pair;
	typedef ListT<Pair> List;

	if (DIMENSION::DIM == 2) {
		assert(pn->data->utp_data !=NULL_PTR);
		List& lpexp = (*CAST(List*, pn->data->utp_data));
		Expression FF[4];
		int countFF[] = { 0, 0, 0, 0 };
		for (typename List::iterator iter = lpexp.begin(); iter != lpexp.end();
				++iter) {
			assert(iter->first->pnode == pn); //
			Face* iterf = iter->first;
			Expression* itere = iter->second;
			int idd = int(iterf->direction) - 4;
			FF[idd].plus(*itere);
			countFF[idd]++;
		}
		for (int i = 0; i < 4; i++) {
			if (countFF[i] > 1) {
				FF[i].times(1.0 / Float(countFF[i])); // get average face value;
			}
		}
		// x direction
		Expression exp_x = FF[2] - FF[0];
		Float l = pn->cell->getD(CSAxis_X);
		Float u = getcVal(pn, u_idx);
		Float CFL_x = u / l;
		//ASSERT_MSG(CFL_x < 0.5, " CFL x > 0.5");
		exp_x.times(-CFL_x);  //negative
		// y direction
		Expression exp_y = FF[1] - FF[3];
		l = pn->cell->getD(CSAxis_Y);
		u = getcVal(pn, v_idx);
		Float CFL_y = u / l;
		//ASSERT_MSG(CFL_y < 0.5, " CFL y > 0.5");
		exp_y.times(-CFL_y); //negative
		//
		exp.plus(exp_x);
		exp.plus(exp_y);

		_clear_utp_data(pn);
	}
	return 1;
}

template<class DIMENSION>
int Advection_Eq<DIMENSION>::advance(int step) {
	typedef typename Forest_::pNode pNode;
	for (int i = 0; i < step; i++) {
		// Prediction step --------------------------------
		//1 Traverse face
		for (typename Forest_::iterator it = pforest->begin();
				it != pforest->end(); ++it) {
			_face_exp_fou(it.get_pointer());
		}
		//2 add time
		for (typename Forest_::iterator it = pforest->begin();
				it != pforest->end(); ++it) {
			pNode pn = it.get_pointer();
			Expression exp;
			_node_exp_adv_t(pn, getcVal(pn, dt_idx) * 0.5, exp); //advance half dt
			exp.Insert(ExpTerm(getIDX(pn), pn, 1.0));
			refcVal(pn, phin_idx) = exp.cal_val(phi_idx);
			//
		}
		// Correction step add limiter: limite variables ------------
		for (typename Forest_::iterator it = pforest->begin();
				it != pforest->end(); ++it) {
			_face_exp_tvd(it.get_pointer(), phin_idx);
		}
		for (typename Forest_::iterator it = pforest->begin();
				it != pforest->end(); ++it) {
			pNode pn = it.get_pointer();
			Expression exp;
			_node_exp_adv_t(pn, getcVal(pn, dt_idx), exp);
			Float valn = exp.cal_val(phin_idx);
			refcVal(pn, phinn_idx) = getcVal(pn, phi_idx) + valn;
		}

		// Put back -------------------------------------------------
		for (typename Forest_::iterator it = pforest->begin();
				it != pforest->end(); ++it) {
			pNode pn = it.get_pointer();
			refcVal(pn, phi_idx) = getcVal(pn, phinn_idx);
		}
	}
	return 0;
}

template<class DIMENSION>
int Advection_Eq<DIMENSION>::_substitute_boudary_val(Expression& exp, st idx) {
	for (typename Expression::iterator iter = exp.begin(); iter != exp.end();
			++iter) {
		if (iter->idx < 0 && iter->idx != ExpTerm::IDX_CONST) {
			vt v = iter->val * getcVal(iter->pnode, idx);
			ExpTerm ct(ExpTerm::IDX_CONST, NULL_PTR, v);
			exp.Insert(ct);
			iter->val = 0.0;
		}
	}
	exp.trim_zero();
	return 1;
}
template<class DIMENSION>
int Advection_Eq<DIMENSION>::_bulid_matrix_fou(MatrixSCR<Float>& mat,
		arrayListV<Float>& b) {
	typedef typename Forest_::pNode pNode;

	//1 Traverse face
	for (typename Forest_::iterator it = pforest->begin(); it != pforest->end();
			++it) {
		_face_exp_fou(it.get_pointer());
	}
	ListT<st> l_rowptr;
	l_rowptr.push_back(0);
	ListT<st> l_colid;
	ListT<Float> l_val;
	ListT<Float> l_b;
	int countnz = 0;
	//2 build matrix
	for (typename Forest_::iterator it = pforest->begin(); it != pforest->end();
			++it) {
		pNode pn = it.get_pointer();
		Expression exp;
		_node_exp_adv(pn, exp);
		_substitute_boudary_val(exp, phi_idx);
		int fconst = 0;
		for (typename Expression::iterator ite = exp.begin(); ite != exp.end();
				++ite) {
			if (ite->idx != ite->IDX_CONST) {
				l_colid.push_back(ite->idx);
				l_val.push_back(ite->val);
				countnz++;
			} else {
				fconst = 1;
				l_b.push_back(-ite->val);    //!!!!! negative added here
			}
		}
		if (fconst == 0) {
			l_b.push_back(0);
		}
		l_rowptr.push_back(countnz);
	}
	//copy list to array  ======================
	st nr = l_rowptr.size() - 1;
	st nz = l_val.size();
	ASSERT(nz == l_colid.size());
	ASSERT(nr <= nz);
	mat.newsize(nr, nr, nz);
	b.reconstruct(l_b.size());
	int i = 0;
	for (typename ListT<st>::iterator it = l_colid.begin(); it != l_colid.end();
			++it) {
		mat.col_ind(i) = (*it);
		i++;
	}
	i = 0;
	for (typename ListT<st>::iterator it = l_rowptr.begin();
			it != l_rowptr.end(); ++it) {
		mat.row_ptr(i) = (*it);
		i++;
	}
	i = 0;
	for (typename ListT<Float>::iterator it = l_val.begin(); it != l_val.end();
			++it) {
		mat.val(i) = (*it);
		i++;
	}
	i = 0;
	for (typename ListT<Float>::iterator it = l_b.begin(); it != l_b.end();
			++it) {
		b[i] = (*it);
		i++;
	}
	return 1;
}

template<class DIMENSION>
int Advection_Eq<DIMENSION>::_bulid_matrix_tvd(MatrixSCR<Float>& mat,
		arrayListV<Float>& b) {
	typedef typename Forest_::pNode pNode;

	//1 Traverse face
	for (typename Forest_::iterator it = pforest->begin(); it != pforest->end();
			++it) {
		_face_exp_tvd(it.get_pointer(), phi_idx);
	}
	ListT<st> l_rowptr;
	l_rowptr.push_back(0);
	ListT<st> l_colid;
	ListT<Float> l_val;
	ListT<Float> l_b;
	int countnz = 0;
	//2 build matrix
	for (typename Forest_::iterator it = pforest->begin(); it != pforest->end();
			++it) {
		pNode pn = it.get_pointer();
		Expression exp;
		_node_exp_adv(pn, exp);
		_substitute_boudary_val(exp, phi_idx);
		int fconst = 0;
		for (typename Expression::iterator ite = exp.begin(); ite != exp.end();
				++ite) {
			if (ite->idx != ite->IDX_CONST) {
				l_colid.push_back(ite->idx);
				l_val.push_back(ite->val);
				countnz++;
			} else {
				fconst = 1;
				l_b.push_back(-ite->val);    //!!!!! negative added here
			}
		}
		if (fconst == 0) {
			l_b.push_back(0);
		}
		l_rowptr.push_back(countnz);
	}
	//copy list to array  ======================
	st nr = l_rowptr.size() - 1;
	st nz = l_val.size();
	ASSERT(nz == l_colid.size());
	ASSERT(nr <= nz);
	mat.newsize(nr, nr, nz);
	b.reconstruct(l_b.size());
	int i = 0;
	for (typename ListT<st>::iterator it = l_colid.begin(); it != l_colid.end();
			++it) {
		mat.col_ind(i) = (*it);
		i++;
	}
	i = 0;
	for (typename ListT<st>::iterator it = l_rowptr.begin();
			it != l_rowptr.end(); ++it) {
		mat.row_ptr(i) = (*it);
		i++;
	}
	i = 0;
	for (typename ListT<Float>::iterator it = l_val.begin(); it != l_val.end();
			++it) {
		mat.val(i) = (*it);
		i++;
	}
	i = 0;
	for (typename ListT<Float>::iterator it = l_b.begin(); it != l_b.end();
			++it) {
		b[i] = (*it);
		i++;
	}
	return 1;
}

template<class DIMENSION>
int Advection_Eq<DIMENSION>::slove(Float tol) {
	MatrixSCR<Float> mat;
	arrayListV<Float> b;
	this->_bulid_matrix_fou(mat, b);
	arrayListV<Float> x(b.size());
	// set up ========
	int max_iter = 1000;
	ListT<Float> lr;	//list residual
	//solver =======================
	int sf = Dia_BiCGSTAB(mat, x, b, max_iter, tol, lr);
	if (sf != 0 && sf != 1) {
		std::cerr << " >! solver failed \n";
		return -1;
	}
	// put the value back
	for (typename Forest_::iterator iter = pforest->begin();
			iter != pforest->end(); ++iter) {
		refcVal(iter, phi_idx) = x[getIDX(iter)];
	}
	//
	for (int ic = 0; ic < 100; ic++) {   // revise check the residual
		MatrixSCR<Float> mat2;
		arrayListV<Float> b2;
		this->_bulid_matrix_tvd(mat2, b2);
		arrayListV<Float> x2(b2.size());
		int max_iter2 = 1000;
		ListT<Float> lr2;
		//solver =======================
		sf = Dia_BiCGSTAB(mat2, x2, b2, max_iter2, tol, lr2);
		if (sf != 0 && sf != 1) {
			std::cerr << " >! solver failed \n";
			return -1;
		}
		// put the value back
		for (typename Forest_::iterator iter = pforest->begin();
				iter != pforest->end(); ++iter) {
			refcVal(iter, phi_idx) = x2[getIDX(iter)];
		}
	}
	return 0;
}
template<class DIMENSION>
int Advection_Eq<DIMENSION>::slove_fou(Float tol) {
	MatrixSCR<Float> mat;
	arrayListV<Float> b;
	this->_bulid_matrix_fou(mat, b);
	arrayListV<Float> x(b.size());
	// set up ========
	int max_iter = 1000;
	ListT<Float> lr;	     //list residual
	//solver =======================
	int sf = Dia_BiCGSTAB(mat, x, b, max_iter, tol, lr);
	if (sf != 0 && sf != 1) {
		std::cerr << " >! solver failed \n";
		return -1;
	}
	// put the value back
	for (typename Forest_::iterator iter = pforest->begin();
			iter != pforest->end(); ++iter) {
		refcVal(iter, phi_idx) = x[getIDX(iter)];
	}
	return 0;
}
template<class DIMENSION>
int Advection_Eq<DIMENSION>::_clear_utp_data(pNode pn) {
	if (pn->data->utp_data != NULL_PTR) {
		typedef typename Forest_::Face Face;
		typedef Pair<Face*, Expression*> Pair;
		typedef ListT<Pair> List;

		List& lpexp = (*CAST(List*, pn->data->utp_data));
		for (typename List::iterator iter = lpexp.begin(); iter != lpexp.end();
				++iter) {
			delete iter->first;
			delete iter->second;
		}
		delete &lpexp;
		pn->data->utp_data = NULL_PTR;
		return 1;
	}
	return 0;
}
template<class DIMENSION>
int Advection_Eq<DIMENSION>::_find_C(pFace pface, pNode& pC, Float veo_f,
		BCM* pBCM) {
	if (isPlus(pface->direction)) {
		if (veo_f > 0) {
			pC = pface->pnode;  //o
		} else {
			pC = pBCM->find_ghost(getIDX(pface->pnode), pface->direction);
		}
	} else {
		if (veo_f > 0) {
			pC = pBCM->find_ghost(getIDX(pface->pnode), pface->direction);
		} else {
			pC = pface->pnode;
		}
	}
	return 1;
}
template<class DIMENSION>
int Advection_Eq<DIMENSION>::_find_C(pFace pface, pNode& pC, Float veo_f) {
	if (isPlus(pface->direction)) {
		if (veo_f > 0) {
			pC = pface->pnode;  //o
		} else {
			pC = pface->pneighbor;
		}
	} else {
		if (veo_f > 0) {
			pC = pface->pneighbor;
		} else {
			pC = pface->pnode;
		}
	}
	return 1;
}

// first order up wind ------------------
template<class DIMENSION>
int Advection_Eq<DIMENSION>::_face_scheme_boundary_fou(pFace pface,
		Expression& exp) {
	int arr_veo_idx[] = { u_idx, v_idx, u_idx, v_idx, w_idx, w_idx };
	st veo_idx = arr_veo_idx[int(pface->direction) - 4];
	Float veo_f = interpolate_1order_on_face((*pface), veo_idx);
	//get U C D ------------------------------
	pNode pC = NULL_PTR;
	_find_C(pface, pC, veo_f, pBCM);
	//
	// exp should be empty before insert
	exp.Insert(ExpTerm(getIDX(pC), pC, 1.0));
	return 1;
}
template<class DIMENSION>
int Advection_Eq<DIMENSION>::_face_scheme_equal_fou(pFace pface,
		Expression& exp) {
	int arr_veo_idx[] = { u_idx, v_idx, u_idx, v_idx, w_idx, w_idx };
	st veo_idx = arr_veo_idx[int(pface->direction) - 4];
	Float veo_f = interpolate_1order_on_face((*pface), veo_idx);
	//get U C D ------------------------------
	pNode pC = NULL_PTR;
	_find_C(pface, pC, veo_f);
	//
	// exp should be empty before insert
	exp.Insert(ExpTerm(getIDX(pC), pC, 1.0));
	return 1;
}
template<class DIMENSION>
int Advection_Eq<DIMENSION>::_face_scheme_fine_coarse_fou(pFace pface,
		Expression& exp) {
	// This function is same with equal ==========================
	int arr_veo_idx[] = { u_idx, v_idx, u_idx, v_idx, w_idx, w_idx };
	st veo_idx = arr_veo_idx[int(pface->direction) - 4];
	Float veo_f = interpolate_1order_on_face((*pface), veo_idx);
	//get U C D ------------------------------
	pNode pC = NULL_PTR;
	_find_C(pface, pC, veo_f);
	//
	// exp should be empty before insert
	exp.Insert(ExpTerm(getIDX(pC), pC, 1.0));
	return 1;
}
template<class DIMENSION>
int Advection_Eq<DIMENSION>::_face_scheme_fou(pFace pface, Expression& exp) {
	//face type
	switch (pface->face_type) {
	case SPFT_Error:
		break;
	case SPFT_Boundary:
		this->_face_scheme_boundary_fou(pface, exp);
		break;
	case SPFT_Equal: {
		this->_face_scheme_equal_fou(pface, exp);
		break;
	}
	case SPFT_FineCoarse:
		this->_face_scheme_fine_coarse_fou(pface, exp);
		break;
	case SPFT_CoarseFine:
		std::cout << "Corase Fine function unfinish first order upwind\n";
		break;
	default:
		return -1;
	}
	return pface->face_type;
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
