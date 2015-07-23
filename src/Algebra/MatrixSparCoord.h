/************************
 //  \file   MatrixSCO.h
 //  \brief
 // 
 //  \author zhou
 //  \date   16 avr. 2014 
 ***********************/
#ifndef MATRIXSPARCOORD_H_
#define MATRIXSPARCOORD_H_

#include "../TypeDef.h"
#include "../Utility/ArrayList.h"
#include "Matrix.h"

namespace Larus {

template<class VALUE> class MatrixSCC;
template<class VALUE> class MatrixSCR;

template<class VALUE>
class MatrixSCO {
public:
	typedef LarusDef::size_type st;
	typedef const LarusDef::size_type const_st;
private:
	arrayListV<VALUE> val_;    // data values (nz_ elements)
	arrayListV<st> rowind_;    // row_ind (nz_ elements)
	arrayListV<st> colind_;    // col_ind (nz_ elements)

	st nz_;                   // number of nonzeros
	st dim_[2];               // number of rows, cols
public:
	MatrixSCO(void);
	MatrixSCO(const Matrix &M);
	MatrixSCO(const MatrixSCO &S);
	MatrixSCO(st M, st N, st nz, VALUE *val, st *r, st *c);
	MatrixSCO(const MatrixSCC<VALUE> &C);
	MatrixSCO(const MatrixSCR<VALUE> &R);

	VALUE& val(st i);
	st& row_ind(st i);
	st& col_ind(st i);
	const VALUE& val(st i) const;
	const st& row_ind(st i) const;
	const st& col_ind(st i) const;

	st size() const;
	st iLen() const;
	st jLen() const;
	st NumNonzeros() const;

	MatrixSCO& operator=(const MatrixSCO &C);
	MatrixSCO& newsize(st M, st N, st nz);
	MatrixSCO& resize(st M, st N, st nz);

	VALUE operator()(st i, st j) const;

	arrayListV<VALUE> operator*(const arrayListV<VALUE> &x) const;
	arrayListV<VALUE> transMult(const arrayListV<VALUE> &x) const;

	VALUE max()const;
	VALUE min()const;

	void show(st a) const;
};


template<class VALUE>
MatrixSCO<VALUE>::MatrixSCO(void) :
		val_(0), rowind_(0), colind_(0), nz_(0) {
	dim_[0] = 0;
	dim_[1] = 0;
}
template<class VALUE>
MatrixSCO<VALUE>::MatrixSCO(const Matrix &M) {
	dim_[0] = M.iLen();
	dim_[1] = M.jLen();
	// count non-zero
	st countnz = 0;
	for (st i = 0; i < M.iLen(); i++) {
		for (st j = 0; j < M.jLen(); j++) {
			if (M[i][j] != 0) {
				countnz++;
			}
		}
	}
	nz_ = countnz;
	val_.reconstruct(nz_);
	rowind_.reconstruct(nz_);
	colind_.reconstruct(nz_);
	st idx = 0;
	for (st i = 0; i < M.iLen(); i++) {
		for (st j = 0; j < M.jLen(); j++) {
			if (M[i][j] != 0) {
				val_[idx] = M[i][j];
				rowind_[idx] = i;
				colind_[idx] = j;
				idx++;
			}
		}
	}
}
template<class VALUE>
MatrixSCO<VALUE>::MatrixSCO(const MatrixSCO &S) :
		val_(S.val_), rowind_(S.rowind_), colind_(S.colind_), nz_(S.nz_) {
	dim_[0] = S.dim_[0];
	dim_[1] = S.dim_[1];
}
template<class VALUE>
MatrixSCO<VALUE>::MatrixSCO(st M, st N, st nz, VALUE *val, st *r,
		st *c) :
		val_(val, nz), rowind_((*r), nz), colind_(*c, nz), nz_(nz) {
	dim_[0] = M;
	dim_[1] = N;
}
template<class VALUE>
MatrixSCO<VALUE>::MatrixSCO(const MatrixSCR<VALUE> &R) :
		val_(R.NumNonzeros()), rowind_(R.NumNonzeros()), colind_(
				R.NumNonzeros()), nz_(R.NumNonzeros()) {
	dim_[0] = R.getiLen();
	dim_[1] = R.getjLen();
	int count = 0;
	int i, j;
//  Loop through rows...
	for (i = 1; i <= dim_[0]; i++) {
		for (j = count; j < R.row_ptr(i); j++) {
			val_[count] = R.val(count);
			colind_[count] = R.col_ind(count);
			rowind_[count] = i - 1;
			count++;
		}
	}
}
template<class VALUE>
MatrixSCO<VALUE>::MatrixSCO(const MatrixSCC<VALUE> &C) :
		val_(C.NumNonzeros()), rowind_(C.NumNonzeros()), colind_(
				C.NumNonzeros()), nz_(C.NumNonzeros()) {
	dim_[0] = C.getiLen();
	dim_[1] = C.getjLen();

	int count = 0;
	int i, j;
//  Loop through columns...
	for (j = 1; j <= dim_[1]; j++) {
		for (i = count; i < C.col_ptr(j); i++) {
			val_[count] = C.val(count);
			rowind_[count] = C.row_ind(count);
			colind_[count] = j - 1;
			count++;
		}
	}
}

template<class VALUE>
MatrixSCO<VALUE>& MatrixSCO<VALUE>::operator=(const MatrixSCO<VALUE> &C) {
	dim_[0] = C.dim_[0];
	dim_[1] = C.dim_[1];
	nz_ = C.nz_;
	val_ = C.val_;
	rowind_ = C.rowind_;
	colind_ = C.colind_;
	return *this;
}
template<class VALUE>
MatrixSCO<VALUE>& MatrixSCO<VALUE>::newsize(MatrixSCO<VALUE>::st M, MatrixSCO<VALUE>::st N, MatrixSCO<VALUE>::st nz) {
	dim_[0] = M;
	dim_[1] = N;
	nz_ = nz;
	val_.reconstruct(nz);
	rowind_.reconstruct(nz);
	colind_.reconstruct(nz);
	return *this;
}

template<class VALUE>
MatrixSCO<VALUE>& MatrixSCO<VALUE>::resize(MatrixSCO<VALUE>::st M, MatrixSCO<VALUE>::st N, MatrixSCO<VALUE>::st nz) {
	dim_[0] = M;
	dim_[1] = N;
	nz_ = nz;
	val_.resize(nz);
	rowind_.resize(nz);
	colind_.resize(nz);
	return *this;
}

//slow---
template<class VALUE>
VALUE MatrixSCO<VALUE>::operator()(MatrixSCO<VALUE>::st i, MatrixSCO<VALUE>::st j) const {
	assert(i >= 0 && i < dim_[0]);
	assert(j >= 0 && j < dim_[1]);
	for (MatrixSCO<VALUE>::st t = 0; t < nz_; t++) {
		if (rowind_(t) == i && colind_(t) == j) {
			return val_(t);
		}
	}
	return 0.0;
}
template<class VALUE>
arrayListV<VALUE> MatrixSCO<VALUE>::operator*(const arrayListV<VALUE> &x) const {
	st M = dim_[0];
	st N = dim_[1];
//  Check for compatible dimensions:
	assert(x.size() == N);
	arrayListV<VALUE> res(M);
	for (st i = 0; i < nz_; i++) {
		res[rowind_[i]] += x[colind_[i]] * val_[i];
	}
	return res;
}
template<class VALUE>
arrayListV<VALUE> MatrixSCO<VALUE>::transMult(const arrayListV<VALUE> &x) const {
	st tM = dim_[1];
	st tN = dim_[0];
	assert(!(x.Len() == tN));
	arrayListV<VALUE> res(tM);
	for (st i = 0; i < nz_; i++) {
		res[colind_[i]] += x[rowind_[i]] * val_[i];
	}
	return res;
}
template<class VALUE>
VALUE MatrixSCO<VALUE>::max()const{
	return val_.findMax();
}
template<class VALUE>
VALUE MatrixSCO<VALUE>::min()const{
	return val_.findMin();
}

template<class VALUE>
void MatrixSCO<VALUE>::show(MatrixSCO<VALUE>::st a) const {
	//if a==0 show matrix in Coordinate formate
	if (a == 0) {
		std::cout << "RowIdx " << "ColIdx " << "Value " << std::endl;
		for (st i = 0; i < nz_; i++) {
			std::cout << std::scientific << rowind_[i] << "  ";
			std::cout << std::scientific << colind_[i] << "  ";
			std::cout << std::scientific << val_[i] << "\n";
		}
	} else {
		for (st i = 0; i < this->iLen(); i++) {
			for (st j = 0; j < this->jLen(); j++) {
				bool isnz = 0;
				for (st t = 0; t < nz_; t++) {
					if (rowind_(t) == i && colind_(t) == j) {
						std::cout << std::scientific << val_(t) << "  ";
						isnz = 1;
					}
				}
				if (isnz == 0) {
					std::cout << std::scientific << 0.0 << "  ";
				}
			}
			std::cout << std::endl;
		}
	}

}
template<class VALUE>
VALUE& MatrixSCO<VALUE>::val(st i) {
	return val_(i);
}
template<class VALUE>
typename MatrixSCO<VALUE>::st& MatrixSCO<VALUE>::row_ind(MatrixSCO<VALUE>::st i) {
	return rowind_(i);
}
template<class VALUE>
typename MatrixSCO<VALUE>::st& MatrixSCO<VALUE>::col_ind(MatrixSCO<VALUE>::st i) {
	return colind_(i);
}
template<class VALUE>
const VALUE& MatrixSCO<VALUE>::val(MatrixSCO<VALUE>::st i) const {
	return val_(i);
}
template<class VALUE>
typename MatrixSCO<VALUE>::const_st& MatrixSCO<VALUE>::row_ind(MatrixSCO<VALUE>::st i) const {
	return rowind_(i);
}
template<class VALUE>
typename MatrixSCO<VALUE>::const_st& MatrixSCO<VALUE>::col_ind(MatrixSCO<VALUE>::st i) const {
	return colind_(i);
}
template<class VALUE>
typename MatrixSCO<VALUE>::st MatrixSCO<VALUE>::size() const {
	return dim_[0] * dim_[1];
}
template<class VALUE>
typename MatrixSCO<VALUE>::st MatrixSCO<VALUE>::iLen() const {
	return dim_[0];
}
template<class VALUE>
typename MatrixSCO<VALUE>::st MatrixSCO<VALUE>::jLen() const {
	return dim_[1];
}
template<class VALUE>
typename MatrixSCO<VALUE>::st MatrixSCO<VALUE>::NumNonzeros() const {
	return nz_;
}



}
// This is the end of namespace

#endif /* MATRIXSPARCOORD_H_ */
