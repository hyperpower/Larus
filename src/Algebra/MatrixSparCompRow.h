/************************
 //  \file   MatrixSCR.h
 //  \brief
 // 
 //  \author zhou
 //  \date   25 avr. 2014 
 ***********************/
#ifndef MATRIXSPARCOMPROW_H_
#define MATRIXSPARCOMPROW_H_

#include "../TypeDef.h"
#include "../Utility/ArrayList.h"
#include "Matrix.h"

namespace Larus {

template<class VALUE> class MatrixSCO;
template<class VALUE> class MatrixSCC;

template<class VALUE>
class MatrixSCR {
public:
	typedef LarusDef::size_type st;
	typedef const LarusDef::size_type const_st;
private:
	arrayListV<VALUE> val_;    // data values (nz_ elements)
	arrayListV<st> rowptr_;    // row_ptr     (dim_[0]+1 elements)
	arrayListV<st> colind_;    // col_ind     (nz_ elements)

	st nz_;                   // number of nonzeros
	st dim_[2];               // number of rows, cols

public:
	MatrixSCR(void);
	//MatrixSCC(const Matrix &M);
	MatrixSCR(const MatrixSCR<VALUE> &S);
	MatrixSCR(const MatrixSCC<VALUE> &C);
	MatrixSCR(const MatrixSCO<VALUE> &CO);
	MatrixSCR(st M, st N, st nz, VALUE *val, st *r, st *c);
	MatrixSCR(st M, st N, st nz, const arrayListV<VALUE> &val,
			const arrayListV<st> &r, const arrayListV<st> &c);

	VALUE& val(st i);
	st& col_ind(st i);
	st& row_ptr(st i);
	const VALUE& val(st i) const;
	const st& col_ind(st i) const;
	const st& row_ptr(st i) const;

	st size() const;
	st iLen() const;
	st jLen() const;
	st NumNonzeros() const;

	MatrixSCR& operator=(const MatrixSCR &C);
	MatrixSCR& newsize(st M, st N, st nz);

	VALUE operator()(st i, st j) const;

	arrayListV<VALUE> operator*(const arrayListV<VALUE> &x) const;
	//Matrix operator*(const MatrixSCR &x) const;
	arrayListV<VALUE> transMult(const arrayListV<VALUE> &x) const;

	void trans();
	// slow
	MatrixSCR getTrans() const;
	//
	VALUE max()const;
	VALUE min()const;
	//


	// IO =======================================
	void show(st a) const;
};

template<class VALUE>
MatrixSCR<VALUE>::MatrixSCR(void) :
		val_(0), rowptr_(0), colind_(0), nz_(0) {
	dim_[0] = 0;
	dim_[1] = 0;
}
template<class VALUE>
MatrixSCR<VALUE>::MatrixSCR(const MatrixSCR<VALUE> &S) :
		val_(S.val_), rowptr_(S.rowptr_), colind_(S.colind_), nz_(S.nz_) {
	dim_[0] = S.dim_[0];
	dim_[1] = S.dim_[1];
}
template<class VALUE>
MatrixSCR<VALUE>::MatrixSCR(const MatrixSCC<VALUE> &C) :
		val_(C.NumNonzeros()),  //
		rowptr_(C.iLen() + 1),  //
		colind_(C.NumNonzeros()),  //
		nz_(C.NumNonzeros())  //
{
	dim_[0] = C.getiLen();
	dim_[1] = C.getjLen();

	st i, j;

	arrayListV<st> tally(C.iLen() + 1, 0);
	//      First pass through nonzeros.  Tally entries in each row.
	//      And calculate rowptr array.
	for (i = 0; i < nz_; i++) {
		tally[C.row_ind(i)]++;
	}
	rowptr_[0] = 0;
	for (j = 0; j < dim_[0]; j++)
		rowptr_(j + 1) = rowptr_(j) + tally(j);
	//      Make copy of rowptr for use in second pass.
	tally = rowptr_;
	//      Second pass through nonzeros.   Fill in index and value entries.
	st count = 0;
	for (i = 1; i <= dim_[1]; i++) {
		for (j = count; j < C.col_ptr(i); j++) {
			val_[tally(C.row_ind(j))] = C.val(j);
			colind_[tally(C.row_ind(j))] = i - 1;
			tally[C.row_ind(count)]++;
			count++;
		}
	}
}
template<class VALUE>
MatrixSCR<VALUE>::MatrixSCR(const MatrixSCO<VALUE> &CO) :
		val_(CO.NumNonzeros()), rowptr_(CO.iLen() + 1), colind_(
				CO.NumNonzeros()), nz_(CO.NumNonzeros()) {
	dim_[0] = CO.iLen();
	dim_[1] = CO.jLen();

	st i;
	arrayListV<st> tally(CO.iLen() + 1, 0);
	//      First pass through nonzeros.  Tally entries in each row.
	//      And calculate rowptr array.
	for (i = 0; i < nz_; i++) {
		tally[CO.row_ind(i)]++;
	}
	rowptr_(0) = 0;
	for (i = 0; i < dim_[0]; i++) {
		rowptr_(i + 1) = rowptr_(i) + tally(i);
	}
	//      Make copy of rowptr for use in second pass.
	tally = rowptr_;
	//      Second pass through nonzeros.   Fill in index and value entries.
	for (i = 0; i < nz_; i++) {
		val_[tally(CO.row_ind(i))] = CO.val(i);
		colind_[tally(CO.row_ind(i))] = CO.col_ind(i);
		tally[CO.row_ind(i)]++;
	}
}
template<class VALUE>
MatrixSCR<VALUE>::MatrixSCR(typename MatrixSCR<VALUE>::st M,   //
		typename MatrixSCR<VALUE>::st N,   //
		typename MatrixSCR<VALUE>::st nz,  //
		VALUE *val,                        //
		typename MatrixSCR<VALUE>::st *r,  //
		typename MatrixSCR<VALUE>::st *c) :
		val_(val, nz), rowptr_(*r, M + 1), colind_(*c, nz), nz_(nz) {
	dim_[0] = M;
	dim_[1] = N;
}
template<class VALUE>
MatrixSCR<VALUE>::MatrixSCR(typename MatrixSCR<VALUE>::st M,  //
		typename MatrixSCR<VALUE>::st N,  //
		typename MatrixSCR<VALUE>::st nz, //
		const arrayListV<VALUE> &val,     //
		const arrayListV<st> &r,          //
		const arrayListV<st> &c) :
		val_(val), rowptr_(r), colind_(c), nz_(nz) {
	dim_[0] = M;
	dim_[1] = N;
}
template<class VALUE>
VALUE& MatrixSCR<VALUE>::val(typename MatrixSCR<VALUE>::st i) {
	return val_(i);
}
template<class VALUE>
typename MatrixSCR<VALUE>::st& MatrixSCR<VALUE>::col_ind(
		typename MatrixSCR<VALUE>::st i) {
	return colind_(i);
}
template<class VALUE>
typename MatrixSCR<VALUE>::st& MatrixSCR<VALUE>::row_ptr(
		typename MatrixSCR<VALUE>::st i) {
	return rowptr_(i);
}
template<class VALUE>
const VALUE& MatrixSCR<VALUE>::val(typename MatrixSCR<VALUE>::st i) const {
	return val_(i);
}
template<class VALUE>
const typename MatrixSCR<VALUE>::st& MatrixSCR<VALUE>::col_ind(
		typename MatrixSCR<VALUE>::st i) const {
	return colind_(i);
}
template<class VALUE>
const typename MatrixSCR<VALUE>::st& MatrixSCR<VALUE>::row_ptr(
		typename MatrixSCR<VALUE>::st i) const {
	return rowptr_(i);
}
template<class VALUE>
typename MatrixSCR<VALUE>::st MatrixSCR<VALUE>::size() const {
	return dim_[0] * dim_[1];
}
template<class VALUE>
typename MatrixSCR<VALUE>::st MatrixSCR<VALUE>::iLen() const {
	return dim_[0];
}
template<class VALUE>
typename MatrixSCR<VALUE>::st MatrixSCR<VALUE>::jLen() const {
	return dim_[1];
}
template<class VALUE>
typename MatrixSCR<VALUE>::st MatrixSCR<VALUE>::NumNonzeros() const {
	return nz_;
}
template<class VALUE>
VALUE MatrixSCR<VALUE>::max()const{
	return val_.findMax();
}
template<class VALUE>
VALUE MatrixSCR<VALUE>::min()const{
	return val_.findMin();
}

template<class VALUE>
MatrixSCR<VALUE>& MatrixSCR<VALUE>::operator=(const MatrixSCR<VALUE> &R) {
	dim_[0] = R.dim_[0];
	dim_[1] = R.dim_[1];
	nz_ = R.nz_;
	val_ = R.val_;
	rowptr_ = R.rowptr_;
	colind_ = R.colind_;
	return *this;
}
template<class VALUE>
MatrixSCR<VALUE>& MatrixSCR<VALUE>::newsize(typename MatrixSCR<VALUE>::st M,
		typename MatrixSCR<VALUE>::st N, typename MatrixSCR<VALUE>::st nz) {
	dim_[0] = M;
	dim_[1] = N;
	nz_ = nz;
	val_.reconstruct(nz);
	rowptr_.reconstruct(M + 1);
	colind_.reconstruct(nz);
	return *this;
}
template<class VALUE>
VALUE MatrixSCR<VALUE>::operator()(typename MatrixSCR<VALUE>::st i,
		typename MatrixSCR<VALUE>::st j) const {
	assert(i >= 0 && i < dim_[0]);
	assert(j >= 0 && j < dim_[1]);
	for (st t = rowptr_(i); t < rowptr_(i + 1); t++) {
		if (colind_(t) == j)
			return val_(t);
	}
	return 0.0;
}
template<class VALUE>
arrayListV<VALUE> MatrixSCR<VALUE>::operator*(
		const arrayListV<VALUE> &x) const {
	st M = dim_[0];
	st N = dim_[1];
	//  Check for compatible dimensions:
	assert(x.size() == N);

	arrayListV<VALUE> res(M);
	for (st i = 0; i < M; ++i) {
		for (st j = rowptr_[i]; j < rowptr_[i + 1]; ++j) {
			res[i] += x[colind_[j]] * val_[j];
		}
	}
	return res;
}

template<class VALUE>
arrayListV<VALUE> MatrixSCR<VALUE>::transMult(
		const arrayListV<VALUE> &x) const {
	st Mt = dim_[1];
	st Nt = dim_[0];
	//  Check for compatible dimensions:
	assert(x.size() == Nt);
	arrayListV<VALUE> res(Mt);
	for (st i = 0; i < Mt; ++i) {
		for (st j = rowptr_[i]; j < rowptr_[i + 1]; ++j) {
			res[i] += x[colind_[j]] * val_[j];
		}
	}
	return res;
}

template<class VALUE>
void MatrixSCR<VALUE>::trans() {
	MatrixSCC<VALUE> C(dim_[1], dim_[0], nz_, val_, colind_, rowptr_);
	dim_[0] = C.getiLen();
	dim_[1] = C.getjLen();

	st i, j;

	arrayListV<st> tally(C.getiLen() + 1, 0);
	//      First pass through nonzeros.  Tally entries in each row.
	//      And calculate rowptr array.
	for (i = 0; i < nz_; i++) {
		tally[C.row_ind(i)]++;
	}
	rowptr_[0] = 0;
	for (j = 0; j < dim_[0]; j++)
		rowptr_(j + 1) = rowptr_(j) + tally(j);
	//      Make copy of rowptr for use in second pass.
	tally = rowptr_;
	//      Second pass through nonzeros.   Fill in index and value entries.
	st count = 0;
	for (i = 1; i <= dim_[1]; i++) {
		for (j = count; j < C.col_ptr(i); j++) {
			val_[tally(C.row_ind(j))] = C.val(j);
			colind_[tally(C.row_ind(j))] = i - 1;
			tally[C.row_ind(count)]++;
			count++;
		}
	}
}
template<class VALUE>
MatrixSCR<VALUE> MatrixSCR<VALUE>::getTrans() const {
	MatrixSCR res(dim_[0], dim_[1], nz_, val_, rowptr_, colind_);
	res.trans();
	return res;
}
template<class VALUE>
void MatrixSCR<VALUE>::show(typename MatrixSCR<VALUE>::st a) const {
	if (a == 0) {
		std::cout << "RowPtr " << "ColInd " << "Value " << std::endl;
		for (st i = 0; i < dim_[0]; i++) {
			for (st ii = rowptr_[i]; ii < rowptr_[i + 1]; ii++) {
				std::cout << std::scientific << i << "  ";
				std::cout << std::scientific << colind_[ii] << "  ";
				std::cout << std::scientific << val_[ii] << "\n";
			}
		}
	} else {
		for (st i = 0; i < dim_[0]; i++) {
			for (st j = 0; j < dim_[1]; j++) {
				bool flag = 0;
				for (st ii = rowptr_[i]; ii < rowptr_[i + 1]; ii++) {
					if (colind_[ii] == j) {
						std::cout << std::scientific << val_[ii] << "  ";
						flag = 1;
					}
				}
				if (flag == 0) {
					std::cout << std::scientific << 0.0 << "  ";
				}
			}
			std::cout << std::endl;
		}
	}
}

}
// This is the end of namespace

#endif /* MATRIXSPARCOMPROW_H_ */
