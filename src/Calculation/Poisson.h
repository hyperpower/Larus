/*
 * Poisson.h
 *
 *  Created on: May 1, 2015
 *      Author: zhou
 */

#ifndef _POISSON_H_
#define _POISSON_H_

#include "../TypeDef.h"
#include "CalDef.h"
#include "../Grid/SPTree.h"
#include "../Grid/SPTreeNode.h"
#include "../Grid/Forest.h"
#include "../Algebra/Arithmetic.h"
#include "../Algebra/Expression.h"

#include "../Algebra/MatrixSparCompRow.h"
#include "../Algebra/Solver_matrix.h"
namespace Larus {
//This file use to solve poisson equation
//
//   ▽•(bata(x,y)  ▽phi(x,y)  ) = f(x,y)     2D --version
//   ▽•(bata(x,y,z)▽phi(x,y,z)) = f(x,y,z)   3D --version

template<class DIMENSION>
class Poisson_Eq: public ObjectBase {
public:
	typedef typename DIMENSION::value_type vt;
	typedef typename DIMENSION::size_type st;
	typedef Forest<DIMENSION> Forest_;
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

	st beta_idx;
	st phi_idx;
	st f_idx;

	Poisson_Eq(Forest_* pf, BCManager<DIMENSION>* pBCM, st, st, st);

	int _f_term(typename Forest_::pNode pn, Expression& exp);

	void _set_val(st idx, pfun_f pfun);
	void set_f(pfun_f pfun);
	void set_beta(pfun_f pfun);
	//
	int _substitute_boudary_val(Expression& exp, st idx);
	int _node_exp(typename Forest_::pNode pn, Expression& exp);
	int _face_scheme_gradient(typename Forest_::pFace pface, Expression& exp);
	int _face_scheme_equal_gradient(typename Forest_::pFace pface,
			Expression& exp);
	int _face_scheme_fine_coarse_gradient_f(typename Forest_::pFace pface,
			Expression& exp);
	int _face_scheme_fine_coarse_gradient_bf(typename Forest_::pFace pface,
			Expression& exp);
	int _face_scheme_boundary_gradient(typename Forest_::pFace pface,
			Expression& exp);
	int _clear_utp_data(typename Forest_::pNode pn);
	int _node_face_exp(typename Forest_::pNode pn);
	int _bulid_matrix(MatrixSCR<Float>& mat, arrayListV<Float>& b);

	//
	int solve(Float tol);

	//IO
	void show() const;
};

template<class DIMENSION>
int Poisson_Eq<DIMENSION>::_f_term(typename Forest_::pNode pn,
		Expression& exp) {
	ExpTerm f_(ExpTerm::IDX_CONST, NULL_PTR, //
			-(pn->data->aCenterData[f_idx] * pn->cell->volume()));  //
	exp.Insert(f_);
	return 1;
}

template<class DIMENSION>
void Poisson_Eq<DIMENSION>::_set_val(st idx, pfun_f pfun) {
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
void Poisson_Eq<DIMENSION>::set_f(typename Poisson_Eq<DIMENSION>::pfun_f pfun) {
	_set_val(f_idx, pfun);
}

template<class DIMENSION>
void Poisson_Eq<DIMENSION>::set_beta(
		typename Poisson_Eq<DIMENSION>::pfun_f pfun) {
	_set_val(beta_idx, pfun);
}

template<class DIMENSION>
int Poisson_Eq<DIMENSION>::_clear_utp_data(typename Forest_::pNode pn) {
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
int Poisson_Eq<DIMENSION>::_node_exp(typename Forest_::pNode pn,
		Expression& exp) {
	typedef typename Forest_::Face Face;
	typedef Pair<Face*, Expression*> Pair;
	typedef ListT<Pair> List;

	if (DIMENSION::DIM == 2) {
		CSAxis arr_dp[] = { CSAxis_Y, CSAxis_X, CSAxis_Y, CSAxis_X };
		//
		assert(pn->data->utp_data !=NULL_PTR);
		List& lpexp = (*CAST(List*, pn->data->utp_data));
		Expression sumCF[4];
		int countCF[] = { 0, 0, 0, 0 };
		for (typename List::iterator iter = lpexp.begin(); iter != lpexp.end();
				++iter) {
			assert(iter->first->pnode == pn); //
			Face* iterf = iter->first;
			Expression* itere = iter->second;
			if (iterf->face_type == SPFT_Equal
					|| iterf->face_type == SPFT_Boundary
					|| iterf->face_type == SPFT_FineCoarse) {
				CSAxis dp = arr_dp[int(iterf->direction) - 4];
				Float a_f = pn->cell->getD(dp);
				itere->times(a_f);
				exp.plus((*itere));
			} else if (iterf->face_type == SPFT_CoarseFine) {
				sumCF[int(iterf->direction) - 4].plus(*itere);
				countCF[int(iterf->direction) - 4]++;
			}
		}
		for (int i = 0; i < 4; i++) {
			if (countCF[i] > 0) {
				sumCF[i].times(1.0 / Float(countCF[i])); // get average face gradient;
				CSAxis dp = arr_dp[i];
				Float a_f = pn->cell->getD(dp);
				sumCF[i].times(a_f);
				exp.plus(sumCF[i]);
			}
		}
		//
		_f_term(pn, exp);
		_clear_utp_data(pn);
	}
	return 1;
}
template<class DIMENSION>
int Poisson_Eq<DIMENSION>::_face_scheme_gradient(typename Forest_::pFace pface,
		Expression& exp) {
	//face type
	switch (pface->face_type) {
	case SPFT_Error:
		break;
	case SPFT_Boundary:
		this->_face_scheme_boundary_gradient(pface, exp);
		break;
	case SPFT_Equal: {
		this->_face_scheme_equal_gradient(pface, exp);
		break;
	}
	case SPFT_FineCoarse:
		this->_face_scheme_fine_coarse_gradient_bf(pface, exp);
		break;
	case SPFT_CoarseFine:
		std::cout << "Corase Fine function unfinish\n";
		break;
	default:
		return -1;
	}
	return pface->face_type;
}
template<class DIMENSION>
int Poisson_Eq<DIMENSION>::_face_scheme_fine_coarse_gradient_f(
		typename Forest_::pFace pface, Expression& exp) {
	// face direction in 4,5,6,7
	//                  4         5         6         7
	CSAxis arr_dd[] = { CSAxis_X, CSAxis_Y, CSAxis_X, CSAxis_Y };
	CSAxis arr_dp[] = { CSAxis_Y, CSAxis_X, CSAxis_Y, CSAxis_X };
	int arr_signd[] = { -1, 1, 1, -1 };
	CSAxis dd = arr_dd[int(pface->direction) - 4];
	CSAxis dp = arr_dp[int(pface->direction) - 4];
	Float beta_f = interpolate_1order_on_face((*pface), beta_idx);

	Float coe;  // coe  = beta_f / ( X_C - X_N )
	coe = beta_f
			/ (pface->pnode->cell->get(dd, eCPL_C)
					- pface->pneighbor->cell->get(dd, eCPL_C));
	ExpTerm phi_C(int(pface->pnode->data->aCenterData[Idx_IDX]), pface->pnode,
			arr_signd[int(pface->direction) - 4] * coe);
	Expression expnei;
	// on perpendicular direction
	Float dis = pface->pnode->cell->get(dp, eCPL_C)
			- pface->pneighbor->cell->get(dp, eCPL_C);
	//interpolate --------------
	interpolate_expression_on_axis_2f(pface->pneighbor, dp, dis, expnei,
			(*pBCM));
	// substitute boundary phi
	for (typename Expression::iterator iter = expnei.begin();
			iter != expnei.end(); ++iter) {
		if (iter->idx < 0 && iter->idx != ExpTerm::IDX_CONST) {
			vt v = iter->val * iter->pnode->data->aCenterData[phi_idx];
			ExpTerm ct(ExpTerm::IDX_CONST, NULL_PTR, v);
			expnei.Insert(ct);
			iter->val = 0;
		}
	}
	expnei.trim_zero();

	expnei.times(-arr_signd[int(pface->direction) - 4] * coe);
	exp.Insert(phi_C);
	exp.plus(expnei);
	return 1;
}
template<class DIMENSION>
int Poisson_Eq<DIMENSION>::_substitute_boudary_val(Expression& exp, st idx) {
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
}
template<class DIMENSION>
int Poisson_Eq<DIMENSION>::_face_scheme_fine_coarse_gradient_bf(
		typename Forest_::pFace pface, Expression& exp) {

	// face direction in 4,5,6,7
	//                  4         5         6         7
	CSAxis arr_dd[] = { CSAxis_X, CSAxis_Y, CSAxis_X, CSAxis_Y };
	CSAxis arr_dp[] = { CSAxis_Y, CSAxis_X, CSAxis_Y, CSAxis_X };
	int arr_signd[] = { -1, 1, 1, -1 };
	CSAxis dd = arr_dd[int(pface->direction) - 4];
	CSAxis dp = arr_dp[int(pface->direction) - 4];
	pNode pc = pface->pnode;
	// find pb
	pNode pb = NULL_PTR;
	T_find_neighbor_2(pc, oppositeDirection(pface->direction), pb);
	if (pb == NULL_PTR) {
		pb = pBCM->find_ghost(getIDX(pc), pface->direction);
	}
	if (pb == NULL_PTR) {
		return _face_scheme_fine_coarse_gradient_f(pface, exp);
	}
	Point point_b;
	Expression exp_b;
	_exp_2node_relation(pc, pb, oppositeDirection(pface->direction), point_b,
			exp_b);
	// pf
	Expression exp_f;
	// on perpendicular direction
	Float dis = pface->pnode->cell->get(dp, eCPL_C)
			- pface->pneighbor->cell->get(dp, eCPL_C);
	//interpolate --------------
	interpolate_expression_on_axis_3fb(pface->pneighbor, dp, dis, exp_f,
			(*pBCM));
	//
	Float x = pc->getPoint(pface->direction)[dd];
	Float x1 = pc->cell->get(dd, eCPL_C);
	Exp2D y1;
	y1.Insert(ExpTerm2D(getIDX(pc), pc, 1.0));
	Float x2 = point_b[dd];
	Float x3 = pface->pneighbor->cell->get(dd, eCPL_C);
	exp = second_order_gradient(x, x1, y1, x2, exp_b, x3, exp_f);

	// substitute boundary phi
	_substitute_boudary_val(exp, phi_idx);

	Float beta_f = interpolate_1order_on_face((*pface), beta_idx);

	exp.times(arr_signd[int(pface->direction) - 4] * beta_f);
	return 1;
}

//this function do not calculate the face area
template<class DIMENSION>
int Poisson_Eq<DIMENSION>::_face_scheme_equal_gradient(
		typename Forest_::pFace pface, Expression& exp) {
	// face direction in 4,5,6,7
	//                  4         5         6         7
	CSAxis arr_dd[] = { CSAxis_X, CSAxis_Y, CSAxis_X, CSAxis_Y };
	//CSAxis arr_dp[] = { CSAxis_Y, CSAxis_X, CSAxis_Y, CSAxis_X };
	int arr_signd[] = { -1, 1, 1, -1 };
	CSAxis dd = arr_dd[int(pface->direction) - 4];
	//CSAxis dp = arr_dp[int(pface->direction) - 4];
	Float beta_f = interpolate_1order_on_face((*pface), beta_idx);
	//Float a_f = pface->pnode->cell->getD(dp)
	//		* ((Forest_::Dim == 2) ? 1 : pface->pnode->cell->getD(dp));
	Float coe;
	coe = beta_f
			/ (pface->pnode->cell->get(dd, eCPL_C)
					- pface->pneighbor->cell->get(dd, eCPL_C));
	ExpTerm phi_C(int(pface->pnode->data->aCenterData[Idx_IDX]), pface->pnode,
			arr_signd[int(pface->direction) - 4] * coe);
	ExpTerm phi_N(int(pface->pneighbor->data->aCenterData[Idx_IDX]),
			pface->pneighbor, -arr_signd[int(pface->direction) - 4] * coe);
	exp.Insert(phi_C);
	exp.Insert(phi_N);
	return 1;
}

template<class DIMENSION>
int Poisson_Eq<DIMENSION>::_face_scheme_boundary_gradient(
		typename Forest_::pFace pface, Expression& exp) {
// face direction in 4,5,6,7
// pface->pneighbor must equal to NULL
	int nid = pface->pnode->data->aCenterData[Idx_IDX];
	pface->pneighbor = pBCM->find_ghost(nid, pface->direction);
	CSAxis arr_dd[] = { CSAxis_X, CSAxis_Y, CSAxis_X, CSAxis_Y };
	//CSAxis arr_dp[] = { CSAxis_Y, CSAxis_X, CSAxis_Y, CSAxis_X };
	int arr_signd[] = { -1, 1, 1, -1 };
	CSAxis dd = arr_dd[int(pface->direction) - 4];
	//CSAxis dp = arr_dp[int(pface->direction) - 4];
	Float beta_f = interpolate_1order_on_face((*pface), beta_idx);
	//Float a_f = pface->pnode->cell->getD(dp)
	//		* ((Forest2D::Dim == 2) ? 1 : pface->pnode->cell->getD(dp));
	Float coe;
	coe = beta_f
			/ (pface->pnode->cell->get(dd, eCPL_C)
					- pface->pneighbor->cell->get(dd, eCPL_C));
	ExpTerm phi_C(int(pface->pnode->data->aCenterData[Idx_IDX]), pface->pnode,
			arr_signd[int(pface->direction) - 4] * coe);
	ExpTerm phi_N(ExpTerm::IDX_CONST, NULL_PTR,
			-arr_signd[int(pface->direction) - 4] * coe
					* pface->pneighbor->data->aCenterData[phi_idx]);
	exp.Insert(phi_C);
	exp.Insert(phi_N);
	return 1;
}

template<class DIMENSION>
int Poisson_Eq<DIMENSION>::_node_face_exp(typename Forest_::pNode pn) {
	typedef typename Forest_::Face Face;
	typedef Pair<Face*, Expression*> Pair;
	typedef ListT<Pair> List;
	typedef typename Forest_::pNode pNode;
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
				_face_scheme_gradient(f, *pexp);
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
				pexpn->times(-1.0);
				List& lpexp = (*CAST(List*, pnei->data->utp_data));
				Pair pairn(fn, pexpn);
				lpexp.push_back(pairn);
			}
		}
	}
	return 1;
}

template<class DIMENSION>
int Poisson_Eq<DIMENSION>::_bulid_matrix( //
		MatrixSCR<Float>& mat,  //
		arrayListV<Float>& b) { //
	//1 Traverse face
	for (typename Forest_::iterator it = pforest->begin(); it != pforest->end();
			++it) {
		_node_face_exp(it.get_pointer());
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
		Expression exp;
		_node_exp(it.get_pointer(), exp);

		exp.trim_zero();
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
		//exp.show();
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
int Poisson_Eq<DIMENSION>::solve(Float tol) {
	MatrixSCR<Float> mat;
	arrayListV<Float> b;
	this->_bulid_matrix(mat, b);
	arrayListV<Float> x(b.size());
	//gnuplot_show(mat);
//set up ========
	int max_iter = 1000;
	ListT<Float> lr;	//list residual
//solver =======================
	int sf = Dia_BiCGSTAB(mat, x, b, max_iter, tol, lr);
	if (sf != 0) {
		std::cerr << " >! Poisson solve failed \n";
		return -1;
	}
	//gnuplot_show_ylog(lr);
//put the value back
	for (typename Forest_::iterator iter = pforest->begin();
			iter != pforest->end(); ++iter) {
		iter->data->aCenterData[phi_idx] = x[iter->data->aCenterData[Idx_IDX]];
	}
	return 1;
}

template<class DIMENSION>
Poisson_Eq<DIMENSION>::Poisson_Eq(Forest_* pf, BCManager<DIMENSION>* pbcm, st b,
		st p, st f) :
		pforest(pf), pBCM(pbcm), beta_idx(b), phi_idx(p), f_idx(f) {
}

template<class DIMENSION>
void Poisson_Eq<DIMENSION>::show() const {
	typedef typename DIMENSION::CellData::value_type vt;
	ListT<vt> lxc, lyc, lxm, lxp, lym, lyp, lval;
	for (typename Forest_::const_iterator iter = pforest->begin();
			iter != pforest->end(); iter++) {
		const typename Forest_::Node* pnode = iter.get_pointer();
		lxc.push_back(pnode->cell->getCenterPoint().x);
		lyc.push_back(pnode->cell->getCenterPoint().y);
		lxm.push_back(pnode->cell->getMM().x);
		lxp.push_back(pnode->cell->getPP().x);
		lym.push_back(pnode->cell->getMM().y);
		lyp.push_back(pnode->cell->getPP().y);
		lval.push_back(pnode->data->aCenterData[phi_idx]);
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
		gp.set_cbrange(-1, 1);
	} else {
		gp.set_cbrange(min, max);
	}

	std::ostringstream ss;
	gp.set_xlabel(ss.str());
	gp.set_equal_ratio();
	gp.plot_7(lxc, lyc, lxm, lxp, lym, lyp, lval, cmdstr);
}

} //end name space

#endif /* CALCULATION_POISSON_H_ */
