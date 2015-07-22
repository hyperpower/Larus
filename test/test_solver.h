// *
// * test_solver.h
// *
// *  Created on: May 3, 2015
// *      Author: zhou
//

#ifndef _TEST_SOLVER_H_
#define _TEST_SOLVER_H_

#include "../src/IO/IO_gnuplot.h"
#include "../src/IO/IO_vtk.h"
#include "../src/IO/mmio.h"
#include "../src/Geometry/Triangle.h"
#include "../src/Geometry/Line.h"
#include "../src/Utility/Array.h"
#include "../src/Geometry/Relation.h"
#include "../src/Geometry/Plane.h"
#include "../src/Grid/SPTreeNode.h"
#include "../src/Grid/SPTree.h"
#include "../src/Grid/Cell.h"
#include "../src/Algebra/Space.h"
#include "../src/Calculation/VOF.h"
#include "../src/Algebra/Arithmetic.h"
#include "../src/Algebra/Cube.h"
#include "../src/TypeDef.h"
#include "../src/Algebra/Interpolation.h"
#include "../src/Algebra/Expression.h"
#include "../src/Algebra/Solver_matrix.h"
#include "../src/Algebra/MatrixSparCoord.h"
#include "../src/Algebra/MatrixSparCompRow.h"
#include "../src/Algebra/MatrixSparCompCol.h"
#include "../src/Calculation/Exp.h"

namespace Larus
{

#ifdef __APPLE__
const string test_solver_DIR = "/Users/zhou/Documents/MATLAB/matrix_case/";
#else
const string test_solver_DIR =
		"/home/czhou/Gerris/Taylor_data/Expansion/exp-re31eo100r1.4L4/";
#endif

void test_gauss_e()
{
	Matrix A(3, 3);
	A[0][0] = 2;
	A[0][1] = 1;
	A[0][2] = -1; //
	A[1][0] = -3;
	A[1][1] = -1;
	A[1][2] = 2; //
	A[2][0] = -2;
	A[2][1] = 1;
	A[2][2] = 2; //

	arrayList b(3);
	arrayList x(3);
	b[0] = 8;
	b[1] = -11;
	b[2] = -3;

	A.show();
	cout << " =====" << endl;
	b.show();

	solver_gaussian_elimination(A, b);
	A.show();
	b.show();
}

void matrix_io()
{
	int ret_code;
	MM_typecode matcode;
	FILE *f;
	int M, N, nz;
	int i, *I, *J;
	double *val;

	string filename = "mat1.txt";
	f = fopen(filename.c_str(), "r");
	if (f != NULL_PTR) {
		cout << "matrix file \n";
	}
	if (mm_read_banner(f, &matcode) != 0) {
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}
	cout << mm_typecode_to_str(matcode) << endl;
	/*  This is how one can screen matrix types if their application */
	/*  only supports a subset of the Matrix Market data types.      */
	//cout << mm_is_sparse(matcode) << endl;
	//cout << mm_is_dense(matcode) << endl;
	if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)) {
		printf("Sorry, this application does not support ");
		printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(1);
	}

	/* find out size of sparse matrix .... */
	if ((ret_code = mm_read_mtx_array_size(f, &M, &N)) != 0)
		exit(1);

	/* reseve memory for matrices */

	I = (int *) malloc(nz * sizeof(int));
	J = (int *) malloc(nz * sizeof(int));
	val = (double *) malloc(nz * sizeof(double));

	/* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
	/*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
	/*   (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

	for (i = 0; i < M * N; i++) {
		fscanf(f, "%lf\n", &val[i]);
		I[i]--; /* adjust from 1-based to 0-based */
		J[i]--;
	}

	if (f != stdin)
		fclose(f);

	/************************/
	/* now write out matrix */
	/************************/

	mm_write_banner(stdout, matcode);
	mm_write_mtx_array_size(stdout, M, N);
	for (i = 0; i < M * N; i++)
		fprintf(stdout, "%10.5f\n", val[i]);
}

void matrix_read_dense()
{
	MM_typecode matcode;
	int M;
	int N;
	double* val;
	int i;

	string filename = "mat1.txt";
	FILE* f = fopen(filename.c_str(), "r");
	if (f != NULL_PTR) {
		cout << "matrix file \n";
	}
	if (mm_read_banner(f, &matcode) != 0) {
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}

	if (!mm_is_matrix(matcode) && !mm_is_dense(matcode) && !mm_is_real(matcode)) {
		printf("Sorry, this function does not support ");
		printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(-1);
	}
	//read dense matrix
	mm_read_mtx_array_size(f, &M, &N);
	val = (double *) malloc(M * N * sizeof(double));

	for (i = 0; i < M * N; i++) {
		fscanf(f, "%lf\n", &val[i]);
	}

	mm_write_banner(stdout, matcode);
	mm_write_mtx_array_size(stdout, M, N);

	/* NOTE: matrix market files use 1-based indices, i.e. first element
	 of a vector has index 1, not 0.  */

	for (i = 0; i < M * N; i++)
		fprintf(stdout, "%10.3f\n", val[i]);

	string fna = "mat1_a.txt";
	Matrix a;
	mm_read_mtx_dense(fna, a);
	//a.show();
	string fnb = "mat1_b.txt";
	arrayList b;
	mm_read_array(fnb, b);
	string fnx = "mat1_x.txt";
	arrayList x;
	mm_read_array(fnx, x);

	solver_gaussian_elimination(a, b);
	arrayList err = x - b;
	//err.show();
	cout << "norm 1   " << nrm1(err) << endl;
	cout << "norm 2   " << nrm2(err) << endl;
	cout << "norm inf " << nrminf(err) << endl;

}

void matrix_write_sparse()
{
	MM_typecode matcode;
	const int nz = 4;
	const int M = 10;
	const int N = 10;
	int I[nz] =
	{ 0, 4, 2, 8 };
	int J[nz] =
	{ 3, 8, 7, 5 };
	double val[nz] =
	{ 1.1, 2.2, 3.2, 4.4 };
	int i;

	mm_initialize_typecode(&matcode);
	mm_set_matrix(&matcode);
	mm_set_sparse(&matcode);
	mm_set_real(&matcode);

	mm_write_banner(stdout, matcode);
	mm_write_mtx_crd_size(stdout, M, N, nz);

	/* NOTE: matrix market files use 1-based indices, i.e. first element
	 of a vector has index 1, not 0.  */

	for (i = 0; i < nz; i++)
		fprintf(stdout, "%d %d %10.3g\n", I[i] + 1, J[i] + 1, val[i]);
}

void test_sp1()
{
	string workdir= test_solver_DIR;
	MatrixSCO<Float> mf;
	string fn_matrix = "685_bus";
	mm_read_mtx_sparse( workdir + fn_matrix +".mtx", mf);
    cout << mf.NumNonzeros()<<endl;
    cout << "matrix information  oo==============\n";
    cout << "i = " << mf.iLen() << "   j = "<< mf.jLen() << endl;
	//mf.show(0);
	MatrixSCR<Float> mfr(mf);
	cout << "matrix information  ==============\n";
	cout << "i = " << mfr.iLen() << "   j = "<< mfr.jLen() << endl;
    
    cout << "matrix information  ==============\n";
    cout << "  " << mf(0,0) << "    "<< mfr(0,0) << endl;
	arrayListV<Float> b(mfr.iLen());
	b.assign(1);
	arrayListV<Float> x(mfr.iLen());
	x.assign(1);
    
    //arrayListV<Float> ax = mf*x;
    //ax.show();
    cout.precision(4);
    cout << mf(20, 22) <<"   "<< mfr(20, 22)<<" \n";
    
	//set up ========
	int max_iter = 1000;
	Float tol = 1e-6;
	ListT<Float> lr;  //list residual
	//solver =======================
	IC_BiCGSTAB(mfr, x, b, max_iter, tol, lr);
	cout<<"max iter = "<<max_iter<<endl;
	cout<<"tol      = "<<tol<<endl;
	gnuplot_show_ylog(lr);
	cout<<"solver jacobi "<<endl;
	x.assign(1);
	lr.clear();
	max_iter = 100000;
	tol = 1e-6;
	Dia_BiCGSTAB(mfr, x, b, max_iter, tol, lr);
	cout<<"max iter = "<<max_iter<<endl;
	cout<<"tol      = "<<tol<<endl;
	gnuplot_show_ylog(lr);
	//
	//gnuplot_show(mfr);
	//int i=0;
	//for(ListT<Float>::iterator it=lr.begin(); it!=lr.end(); it++){
	//	cout<<i<< "  "<<(*it)/nrm2(b)<<endl;
	//	i++;
	//}

	//cout << max_iter << "  " << tol << endl;
	//to row comp

}

} //end namespace

#endif /* TEST_SOLVER_H_ */
