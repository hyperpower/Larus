/************************
 //  \file   MatrixSCC.h
 //  \brief
 // 
 //  \author zhou
 //  \date   17 avr. 2014 
 ***********************/
#ifndef MATRIXSPARCOMPCOL_H_
#define MATRIXSPARCOMPCOL_H_

#include "../TypeDef.h"
#include "../Utility/ArrayList.h"
#include "Matrix.h"

namespace Larus
{

template<class VALUE> class MatrixSCO;
template<class VALUE> class MatrixSCC;

template<class VALUE>
class MatrixSCC
{
public:
	typedef LarusDef::size_type st;
	typedef const LarusDef::size_type const_st;
private:
	arrayListV<VALUE> val_;             // data values (nz_ elements)
	arrayListV<st> rowind_;    // row_ind     (nz_ elements)
	arrayListV<st> colptr_;    // col_ptr     (dim_[1]+1 elements)

	st nz_;                   // number of nonzeros
	st dim_[2];               // number of rows, cols

public:
	MatrixSCC(void);
	//MatrixSCC(const Matrix &M);
	MatrixSCC(const MatrixSCC<VALUE> &S);
	MatrixSCC(const MatrixSCR<VALUE> &R);
	MatrixSCC(const MatrixSCO<VALUE> &CO);
	MatrixSCC(st M, st N, st nz, VALUE*val, st *r, st *c);
	MatrixSCC(st M, st N, st nz, const arrayListV<VALUE> &val,
			const arrayListV<st> &r, const arrayListV<st> &c);

	VALUE& val(st i);
	st& row_ind(st i);
	st& col_ptr(st i);
	const VALUE& val(st i) const;
	const st& row_ind(st i) const;
	const st& col_ptr(st i) const;

	st size() const;
	st iLen() const;
	st jLen() const;
	st NumNonzeros() const;

	MatrixSCC& operator=(const MatrixSCC &C);
	MatrixSCC& newsize(st M, st N, st nz);

	VALUE operator()(st i, st j) const;

	arrayListV<VALUE> operator*(const arrayListV<VALUE> &x) const;
	arrayListV<VALUE> transMult(const arrayListV<VALUE> &x) const;
	void show(st a) const;
};

template<class VALUE>
MatrixSCC<VALUE>::MatrixSCC(void) :
		val_(0), rowind_(0), colptr_(0), nz_(0)
{
	dim_[0] = 0;
	dim_[1] = 0;
}
template<class VALUE>
MatrixSCC<VALUE>::MatrixSCC(const MatrixSCC<VALUE> &S) :
		val_(S.val_), rowind_(S.rowind_), colptr_(S.colptr_), nz_(S.nz_)
{
	dim_[0] = S.dim_[0];
	dim_[1] = S.dim_[1];
}
template<class VALUE>
MatrixSCC<VALUE>::MatrixSCC(
		typename MatrixSCC<VALUE>::st M,
		typename MatrixSCC<VALUE>::st N,
		typename MatrixSCC<VALUE>::st nz,
		VALUE *val,
		typename MatrixSCC<VALUE>::st *r,
		typename MatrixSCC<VALUE>::st *c) :
		val_(val, nz), rowind_(*r, nz), colptr_(*c, N + 1), nz_(nz)
{
	dim_[0] = M;
	dim_[1] = N;
}
template<class VALUE>
MatrixSCC<VALUE>::MatrixSCC(
		typename MatrixSCC<VALUE>::st M, //
		typename MatrixSCC<VALUE>::st N, //
		typename MatrixSCC<VALUE>::st nz, //
		const arrayListV<VALUE> &val,
		const arrayListV<MatrixSCC<VALUE>::st> &rid,
		const arrayListV<MatrixSCC<VALUE>::st> &cp) :
		val_(val), rowind_(rid), colptr_(cp), nz_(nz)
{
	dim_[0] = M;
	dim_[1] = N;
}
template<class VALUE>
MatrixSCC<VALUE>::MatrixSCC(const MatrixSCR<VALUE> &R) :
		val_(R.NumNonzeros()), rowind_(R.NumNonzeros()), colptr_(
				R.getjLen() + 1), nz_(R.NumNonzeros())
{
	dim_[0] = R.getiLen();
	dim_[1] = R.getjLen();

	int i, j;
	arrayListV<int> tally(R.getjLen() + 1, 0);
	//      First pass through nonzeros.  Tally entries in each column.
	//      And calculate colptr array.
	for (i = 0; i < nz_; i++) {
		tally[R.col_ind(i)]++;
	}
	colptr_[0] = 0;
	for (j = 0; j < dim_[1]; j++)
		colptr_(j + 1) = colptr_(j) + tally(j);
	//      Make copy of colptr for use in second pass.
	tally = colptr_;
	//      Second pass through nonzeros.  Fill in index and value entries.
	int count = 0;
	for (i = 1; i <= dim_[0]; i++) {
		for (j = count; j < R.row_ptr(i); j++) {
			val_[tally(R.col_ind(j))] = R.val(j);
			rowind_[tally(R.col_ind(j))] = i - 1;
			tally[R.col_ind(count)]++;
			count++;
		}
	}
}
template<class VALUE>
MatrixSCC<VALUE>::MatrixSCC(const MatrixSCO<VALUE> &CO) :
		val_(CO.NumNonzeros()), rowind_(CO.NumNonzeros()), colptr_(
				CO.getjLen() + 1), nz_(CO.NumNonzeros())
{

	dim_[0] = CO.getiLen();
	dim_[1] = CO.getjLen();

	int i, j;
	arrayListV<int> tally(CO.getjLen() + 1, 0);
//  First pass through nonzeros.  Tally entries in each column.
//  And calculate colptr array.
	for (i = 0; i < nz_; i++) {
		tally[CO.col_ind(i)]++;
	}
	colptr_[0] = 0;
	for (j = 0; j < dim_[1]; j++) {
		colptr_[j + 1] = colptr_[j] + tally[j];
	}
//  Make copy of colptr for use in second pass.
	tally = colptr_;
//  Second pass through nonzeros.   Fill in index and value entries.
	for (i = 0; i < nz_; i++) {
		val_[tally[CO.col_ind(i)]] = CO.val(i);
		rowind_[tally[CO.col_ind(i)]] = CO.row_ind(i);
		tally[CO.col_ind(i)]++;
	}
}
template<class VALUE>
VALUE& MatrixSCC<VALUE>::val(typename MatrixSCC<VALUE>::st i)
{
	return val_(i);
}
template<class VALUE>
typename MatrixSCC<VALUE>::st& MatrixSCC<VALUE>::row_ind(
		typename MatrixSCC<VALUE>::st i)
{
	return rowind_(i);
}
template<class VALUE>
typename MatrixSCC<VALUE>::st& MatrixSCC<VALUE>::col_ptr(
		typename MatrixSCC<VALUE>::st i)
{
	return colptr_(i);
}
template<class VALUE>
const VALUE& MatrixSCC<VALUE>::val(typename MatrixSCC<VALUE>::st i) const
{
	return val_(i);
}
template<class VALUE>
const typename MatrixSCC<VALUE>::st& MatrixSCC<VALUE>::row_ind(
		typename MatrixSCC<VALUE>::st i) const
{
	return rowind_(i);
}
template<class VALUE>
const typename MatrixSCC<VALUE>::st& MatrixSCC<VALUE>::col_ptr(
		typename MatrixSCC<VALUE>::st i) const
{
	return colptr_(i);
}
template<class VALUE>
typename MatrixSCC<VALUE>::st MatrixSCC<VALUE>::size() const
{
	return dim_[0] * dim_[1];
}
template<class VALUE>
typename MatrixSCC<VALUE>::st MatrixSCC<VALUE>::iLen() const
{
	return dim_[0];
}
template<class VALUE>
typename MatrixSCC<VALUE>::st MatrixSCC<VALUE>::jLen() const
{
	return dim_[1];
}
template<class VALUE>
typename MatrixSCC<VALUE>::st MatrixSCC<VALUE>::NumNonzeros() const
{
	return nz_;
}
template<class VALUE>
MatrixSCC<VALUE>& MatrixSCC<VALUE>::operator=(const MatrixSCC &C)
{
	dim_[0] = C.dim_[0];
	dim_[1] = C.dim_[1];
	nz_ = C.nz_;
	val_ = C.val_;
	rowind_ = C.rowind_;
	colptr_ = C.colptr_;
	return *this;
}
template<class VALUE>
MatrixSCC<VALUE>& MatrixSCC<VALUE>::newsize(typename MatrixSCC<VALUE>::st M,
		typename MatrixSCC<VALUE>::st N, typename MatrixSCC<VALUE>::st nz)
{
	dim_[0] = M;
	dim_[1] = N;

	nz_ = nz;
	val_.reconstruct(nz);
	rowind_.reconstruct(nz);
	colptr_.reconstruct(N + 1);
	return *this;
}
template<class VALUE>
VALUE MatrixSCC<VALUE>::operator()(typename MatrixSCC<VALUE>::st i,
		typename MatrixSCC<VALUE>::st j) const
{
	assert(i >= 0 && i < dim_[0]);
	assert(j >= 0 && j < dim_[1]);
	for (int t = colptr_(j); t < colptr_(j + 1); t++) {
		if (rowind_(t) == i) {
			return val_(t);
		}
	}
	return 0.0;
}
template<class VALUE>
arrayListV<VALUE> MatrixSCC<VALUE>::operator*(const arrayListV<VALUE> &x) const
{
	int M = dim_[0];
	int N = dim_[1];
	//  Check for compatible dimensions:
	assert(x.Len() == N);
	arrayListV<VALUE> res(M);
	for (int i = 0; i < N; i++) {
		for (int j = colptr_[i]; j < colptr_[i + 1]; j++) {
			res[rowind_[j]] += x[i] * val_[j];
		}
	}
	return res;
}
template<class VALUE>
arrayListV<VALUE> MatrixSCC<VALUE>::transMult(const arrayListV<VALUE> &x) const
{
	int Mt = dim_[1];
	int Nt = dim_[0];
	//  Check for compatible dimensions:
	assert(x.Len() == Nt);
	arrayListT<VALUE> res(Mt);
	for (int i = 0; i < Nt; i++) {
		for (int j = colptr_[i]; j < colptr_[i + 1]; j++) {
			res[rowind_[j]] += x[i] * val_[j];
		}
	}
	return res;
}
template<class VALUE>
void MatrixSCC<VALUE>::show(typename MatrixSCC<VALUE>::st a) const
{
	if (a == 0) {
		std::cout << "RowIdx " << "ColPtr " << "Value " << std::endl;
		for (int i = 0; i < dim_[1]; i++) {
			for (int ii = colptr_[i]; ii < colptr_[i + 1]; ii++) {
				std::cout << std::scientific << rowind_[ii] << "  ";
				std::cout << std::scientific << i << "  ";
				std::cout << std::scientific << val_[ii] << "\n";
			}
		}
	} else {
		for (int i = 0; i < dim_[0]; i++) {
			for (int j = 0; j < dim_[1]; j++) {
				bool flag = 0;
				for (int ii = colptr_[j]; ii < colptr_[j + 1]; ii++) {
					if (rowind_[ii] == i) {
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

#endif /* MATRIXSPARCOMPCOL_H_ */
