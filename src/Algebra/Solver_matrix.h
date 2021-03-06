/*
 * Solver_matrix.h
 *
 *  Created on: May 3, 2015
 *      Author: zhou
 */

#ifndef _SOLVER_MATRIX_H_
#define _SOLVER_MATRIX_H_

#include <iostream>
#include <assert.h>
#include "../TypeDef.h"
#include "../Utility/ArrayList.h"
#include "MatrixSparCompRow.h"
#include "BLAS_Level1.h"
#include "Preconditioner.h"

namespace Larus {
/**
 * slove
 * a1 x + b1 y = c1
 * a2 x + b2 y = c2
 *
 */

template<class VALUE>
int solve(VALUE a1, VALUE b1, VALUE c1, VALUE a2, VALUE b2, VALUE c2, VALUE& x,
		VALUE& y) {
	assert(!(a2 == 0 && b2 == 0));
	assert(!(a1 == 0 && b1 == 0));
	if (a2 == 0) {
		y = double(c2) / double(b2);
	} else {
		y = double(c1 - a1 * c2 / a2) / double(b1 - a1 * b2 / a2 + SMALL);
	}
	x = double(c1 - b1 * y) / double(a1 + SMALL);
}
template<class VALUE>
int solver_gaussian_elimination(  //solver
		MatrixV<VALUE>& A,      //the Matrix     [in]  solver will change matrix
		arrayListV<VALUE>& b    //the b          [in][out]  solver will change b
		) {
	//Assert
	ASSERT(A.iLen() == A.jLen());   //the matrix must be n x n;
	ASSERT(b.size() == A.iLen());     //the size b=n
	int n = b.size();
	int i, j, max_i;
	VALUE max, dum;
	// for each variable find pivot row and perform forward substitution
	for (i = 0; i < (n - 1); ++i) {
		//  find the pivot row
		max_i = i;
		max = ABS(A[i][i]);
		for (int ii = i + 1; ii < n; ++ii)
			if ((dum = ABS(A[ii][i])) > max) {
				max = dum;
				max_i = ii;
			}
		if (max == 0.0) {
			return -1;                // the matrix A is singular
		}
		// and if it differs from the current row, interchange the two rows.
		if (max_i != i) {
			for (j = 0; j < n; j++) {
				dum = A[i][j];
				A[i][j] = A[max_i][j];
				A[max_i][j] = dum;
			}
			dum = b[i];
			b[i] = b[max_i];
			b[max_i] = dum;
		}
		// Perform forward substitution
		for (int ii = i + 1; ii < n; ++ii) {
			dum = -A[ii][i] / A[i][i];
			for (j = i + 1; j < n; ++j) {
				A[ii][j] += dum * A[i][j];
			}
			b[ii] += dum * b[i];
		}
	}
	// Perform backward substitution
	for (i = n - 1; i >= 0; i--) {
		if (A[i][i] == 0.0)
			return -1;           // matrix is singular
		dum = b[i];
		for (j = i + 1; j < n; j++) {
			dum -= A[i][j] * b[j];
		}
		b[i] = dum / A[i][i];
	}
	return 0;
}

//spares ==================================================
//=========================================================

// Jacobi solver
template<class VALUE>
int Jacobi(const MatrixSCR<VALUE> &A,    // A  The matrix
		arrayListV<VALUE> &x,        // x
		const arrayListV<VALUE>& b,  // b
		int &max_iter, Float &tol,   // Tolerance
		ListT<Float>& lresid)        // list residual
		{
	Float resid;
	typename arrayListV<VALUE>::size_type n = b.size();
	MatrixSCR<VALUE> T(A);
	arrayListV<VALUE> C(n, 1.0 / SMALL);
	arrayListV<VALUE> newx(n);

	Float normb = nrm2(b);
	arrayListV<VALUE> r = b - A * x;
	//
	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	// construct T
	typedef typename MatrixSCR<VALUE>::st st;
	st M = T.iLen();
	//st N = T.jLen();

	for (st i = 0; i < M; ++i) {
		// find D
		VALUE dv;
		int flag = -1;
		for (st j = T.row_ptr(i); j < T.row_ptr(i + 1); ++j) {
			if (T.col_ind(j) == i) {
				flag = j;
				dv = T.val(j);
				T.val(j) = 0;
				break;
			}
		}
		for (st j = T.row_ptr(i); j < T.row_ptr(i + 1); ++j) {
			if (j == flag) {
				C(i) = b(T.col_ind(j)) / dv;
			} else {
				T.val(j) = -T.val(j) / dv;
			}
		}
	}

	for (int i = 1; i <= max_iter; ++i) {
		//----
		newx = T * x + C;
		r = b - A * newx;
		//----
		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid <= tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}
		x = newx;
	}

	tol = resid;
	return 1;
}

// CG solves the symmetric positive definite linear
// system Ax=b using the Conjugate Gradient method.
template<class VALUE>
int CG(const MatrixSCR<VALUE> &A,    // A  The matrix
		arrayListV<VALUE> &x,        // x
		const arrayListV<VALUE>& b,  // b
		int &max_iter, Float &tol,   // Tolerance
		ListT<Float>& lresid)        // list residual
		{
	Float resid;
	arrayListV<VALUE> p(b.size()), z(b.size()), q(b.size());
	//Array alpha(1), beta(1), rho(1), rho_1(1);
	VALUE alpha, beta, rho, rho_1;

	Float normb = nrm2(b);
	arrayListV<VALUE> r = b - A * x;
	//
	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; ++i) {
		z = r;           //z = M.solve(r);   ===!!!!!
		rho = dot(r, z); //rho(0) = dot(r, z);

		if (i == 1)
			p = z;
		else {
			beta = rho / rho_1; //beta(0) = rho(0) / rho_1(0);
			p = z + beta * p;   //p = z + beta(0) * p;
		}

		q = A * p;
		alpha = rho / dot(p, q);  //alpha(0) = rho(0) / dot(p, q);
		x += alpha * p;           //x += alpha(0) * p;
		r -= alpha * q;           //r -= alpha(0) * q;

		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid <= tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}
		rho_1 = rho;            //rho_1(0) = rho(0);
	}

	tol = resid;
	return 1;
}

template<class VALUE>
int IC_CG(const MatrixSCR<VALUE> &A,  //
		arrayListV<VALUE> &x,         //
		const arrayListV<VALUE>& b,   //
		int &max_iter, Float &tol,    //
		ListT<Float>& lresid)         //
		{
	Float resid;
	arrayListV<VALUE> p, z, q;
	Float alpha, beta, rho, rho_1;
	ICPre<VALUE> preA(A);
	Float normb = nrm2(b);
	arrayListV<VALUE> r = b - A * x;

	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; i++) {
		//if(i<=max_iter){
		z = preA.solve(r);
		//}else{
		//	z = r;                    //z = M.solve(r);   ===!!!!!
		//}
		rho = dot(r, z);             //rho(0) = dot(r, z);

		if (i == 1)
			p = z;
		else {
			beta = rho / rho_1; //beta(0) = rho(0) / rho_1(0);
			p = z + beta * p;   //p = z + beta(0) * p;
		}

		q = A * p;
		alpha = rho / dot(p, q);  //alpha(0) = rho(0) / dot(p, q);
		x += alpha * p;           //x += alpha(0) * p;
		r -= alpha * q;           //r -= alpha(0) * q;

		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid <= tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}

		rho_1 = rho;            //rho_1(0) = rho(0);
	}

	tol = resid;
	return 1;
}

// CGS solves the unsymmetric linear system Ax = b
// using the Conjugate Gradient Squared method
template<class VALUE>
int IC_CGS(const MatrixSCR<VALUE> &A,    // A  The matrix
		arrayListV<VALUE> &x,        // x
		const arrayListV<VALUE>& b,  // b
		int &max_iter,   //max_iter
		Float &tol,   // Tolerance
		ListT<Float>& lresid) {
	ICPre<VALUE> preA(A);
	Float resid;
	VALUE rho_1, rho_2, alpha, beta;
	arrayListV<VALUE> p, phat, q, qhat, vhat, u, uhat;

	Float normb = nrm2(b);
	arrayListV<VALUE> r = b - A * x;
	arrayListV<VALUE> rtilde = r;

	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; i++) {
		rho_1 = dot(rtilde, r);
		if (rho_1 == 0) {
			tol = nrm2(r) / normb;
			return 2;
		}
		if (i == 1) {
			u = r;
			p = u;
		} else {
			beta = rho_1 / rho_2;
			u = r + beta * q;
			p = u + beta * (q + beta * p);
		}
		phat = preA.solve(p);
		vhat = A * phat;
		alpha = rho_1 / dot(rtilde, vhat);
		q = u - alpha * vhat;
		uhat = preA.solve(u + q);
		x += alpha * uhat;
		qhat = A * uhat;
		r -= alpha * qhat;
		rho_2 = rho_1;

		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid < tol) {
			tol = resid;
			max_iter = i;
			return 0; //converged to the desired tolerance tol within maxit iterations.
		}
	}

	tol = resid;
	return 1;         // iterated maxit times but did not converge.
}

//============================================================
// BiCG solves the unsymmetric linear system Ax = b
// using the Preconditioned BiConjugate Gradient method
template<class VALUE>
int IC_BiCG(const MatrixSCR<VALUE> &A,    // A  The matrix
		arrayListV<VALUE> &x,        // x
		const arrayListV<VALUE>& b,  // b
		int &max_iter,   //max_iter
		Float &tol,   // Tolerance
		ListT<Float>& lresid) {
	ICPre<VALUE> preA(A);
	Float resid;
	VALUE rho_1, rho_2, alpha, beta;
	arrayListV<VALUE> z, ztilde, p, ptilde, q, qtilde;

	Float normb = nrm2(b);
	arrayListV<VALUE> r = b - A * x;
	arrayListV<VALUE> rtilde = r;

	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; i++) {
		z = preA.solve(r);
		ztilde = preA.trans_solve(rtilde);
		rho_1 = dot(z, rtilde);
		if (rho_1 == 0) {
			tol = nrm2(r) / normb;
			max_iter = i;
			return 2;
		}
		if (i == 1) {
			p = z;
			ptilde = ztilde;
		} else {
			beta = rho_1 / rho_2;
			p = z + beta * p;
			ptilde = ztilde + beta * ptilde;
		}
		q = A * p;
		qtilde = A.transMult(ptilde);
		alpha = rho_1 / dot(ptilde, q);
		x += alpha * p;
		r -= alpha * q;
		rtilde -= alpha * qtilde;

		rho_2 = rho_1;

		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid < tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}
	}

	tol = resid;
	return 1;
}

//===============================================================
// BiCGSTAB solves the unsymmetric linear system Ax = b
// using the Preconditioned BiConjugate Gradient Stabilized method
template<class VALUE>
int IC_BiCGSTAB( //
		const MatrixSCR<VALUE> &A,    // A  The matrix
		arrayListV<VALUE> &x,        // x
		const arrayListV<VALUE>& b,  // b
		int &max_iter,   //max_iter
		Float &tol,      // Tolerance
		ListT<Float>& lresid) {
	ICPre<VALUE> preA(A);
	Float resid;
	VALUE rho_1, rho_2, alpha, beta, omega;
	arrayListV<VALUE> p, phat, s, shat, t, v;

	Float normb = nrm2(b);
	arrayListV<VALUE> r = b - A * x;
	arrayListV<VALUE> rtilde = r;

	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; i++) {
		rho_1 = dot(rtilde, r);
		if (rho_1 == 0) {
			tol = nrm2(r) / normb;
			return 2;
		}
		if (i == 1)
			p = r;
		else {
			beta = (rho_1 / rho_2) * (alpha / omega);
			p = r + beta * (p - omega * v);
		}
		phat = preA.solve(p);
		v = A * phat;
		alpha = rho_1 / dot(rtilde, v);
		s = r - alpha * v;
		if ((resid = nrm2(s) / normb) < tol) {
			x += alpha * phat;
			tol = resid;
			return 0;
		}
		shat = preA.solve(s);
		t = A * shat;
		omega = dot(t, s) / dot(t, t);
		x += alpha * phat + omega * shat;
		r = s - omega * t;

		rho_2 = rho_1;

		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid < tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}

		if (omega == 0) {
			tol = nrm2(r) / normb;
			return 3;
		}
	}

	tol = resid;
	return 1;
}
template<class VALUE>
int Dia_BiCGSTAB( //
		const MatrixSCR<VALUE> &A,    // A  The matrix
		arrayListV<VALUE> &x,        // x
		const arrayListV<VALUE>& b,  // b
		int &max_iter,   //max_iter
		Float &tol,      // Tolerance
		ListT<Float>& lresid) {
	DiaPre<VALUE> preA(A, 1);
	Float resid;
	VALUE rho_1, rho_2, alpha, beta, omega;
	arrayListV<VALUE> p, phat, s, shat, t, v;

	Float normb = nrm2(b);
	arrayListV<VALUE> r = b - A * x;
	arrayListV<VALUE> rtilde = r;

	if (normb == 0.0)
		normb = 1;

	if ((resid = nrm2(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; i++) {
		rho_1 = dot(rtilde, r);          //dot(v,v)
		if (rho_1 == 0) {
			tol = nrm2(r) / normb;       //norm(v) / s
			return 2;
		}
		if (i == 1)
			p = r;                       // v=v
		else {
			beta = (rho_1 / rho_2) * (alpha / omega);   // s
			p = r + beta * (p - omega * v);             // v + s*(v-s*v)
		}
		phat = preA.solve(p);
		v = A * phat;
		alpha = rho_1 / dot(rtilde, v);
		s = r - alpha * v;
		if ((resid = nrm2(s) / normb) < tol) {
			x += alpha * phat;
			tol = resid;
			return 0;
		}
		shat = preA.solve(s);
		t = A * shat;
		omega = dot(t, s) / dot(t, t);
		x += alpha * phat + omega * shat;
		r = s - omega * t;

		rho_2 = rho_1;

		resid = nrm2(r) / normb;
		lresid.push_back(resid);
		if (resid < tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}

		if (omega == 0) {
			tol = nrm2(r) / normb;
			return 3;
		}
	}

	tol = resid;
	return 1;
}
} //end namespace

#endif /* ALGEBRA_SOLVER_MATRIX_H_ */
