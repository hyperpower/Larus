/*
 * test_python.h
 *
 *  Created on: May 13, 2015
 *      Author: zhou
 */

#ifndef TEST_PYTHON_H_
#define TEST_PYTHON_H_

#include <iostream>
#include <Python.h>

using namespace std;

namespace Larus {

void main_2D(PyObject* pParse);

int test_python() {
	PyObject* pModule = NULL;
	PyObject* pFunc = NULL;
	PyObject* pParam = NULL;
	PyObject* pResult = NULL;
	const char* pBuffer = NULL;
	int iBufferSize = 0;

	Py_Initialize();        //Initialize
	char *path = (char *) "./";
	PySys_SetPath(path);
	pModule = PyImport_ImportModule("test_python");

	if (!pModule) {
		cout << "get module failed!" << endl;
		exit(0);
	} else {
		cout << "get module success!" << endl;
	}

	string attr_name = "var";
	pFunc = PyObject_GetAttrString(pModule, attr_name.c_str());
	if (!pFunc) {
		cout << "get func failed!" << endl;
		//cout << pFunc << endl;
		exit(0);
	} else {
		cout << PyList_Size(pFunc) << endl;
		PyObject* t0 = PyList_GetItem(pFunc, 0);
		cout << t0 << endl;
		cout << PyString_Check(t0) << endl;
		cout << "get func success!" << endl;
	}
	pParam = Py_BuildValue("(s)", "HEHEHE");
	pResult = PyEval_CallObject(pFunc, pParam);

	if (pResult) {
		if (PyArg_Parse(pResult, "(si)", &pBuffer, &iBufferSize)) {
			cout << "=============!" << endl;
			//cout << pBuffer << endl;
			cout << iBufferSize << endl;
		} else {
			cout << "PyArg_Parse erro!" << endl;
		}
	}
	Py_DECREF(pParam);
	Py_DECREF(pFunc);

	Py_Finalize();
	//cout << "hello" << endl;
	return 0;
}

int get_complete_dim(PyObject* model) {
	PyObject* pdim = NULL;
	string attr_name = "_complete_cs_";
	pdim = PyObject_GetAttrString(model, attr_name.c_str());
	if (!pdim) {
		cerr << " >! Cann't find _complete_cs_" << endl;
		exit(0);
	}
	if (!PyInt_Check(pdim)) {
		cerr << " >! Cann't find dim type error" << endl;
		exit(0);
	}
	return int(PyFloat_AsDouble(pdim));
}

inline int _pytype_to_int(PyObject* t, const string& err_msg) {
	if (!PyInt_Check(t)) {
		cerr << err_msg << " type=" << t->ob_type->tp_name << " should be int"
				<< endl;
		exit(0);
	}
	return int(PyFloat_AsDouble(t));
}
inline Float _pytype_to_float(PyObject* t, const string& err_msg) {
	if (!PyFloat_Check(t)) {
		cerr << err_msg << " type=" << t->ob_type->tp_name << " should be float"
				<< endl;
		exit(0);
	}
	return Float(PyFloat_AsDouble(t));
}
inline string _pytype_to_string(PyObject* t, const string& err_msg) {
	if (1 != PyString_Check(t)) {
		cerr << err_msg << endl;
		exit(0);
	}
	return PyString_AsString(t);
}
//the parameter must be as same as forest
void py_construct_param_forest(
		PyObject* model,  //
		int& i, int& j, Float& ox, Float& oy, Float& length, int& maxl, int& k,
		Float& oz) {
	PyObject* pforest = NULL;
	string attr_name = "_complete_forest_";
	pforest = PyObject_GetAttrString(model, attr_name.c_str());
	if (!pforest) {
		cerr << " >! Can't find _complete_forest_" << endl;
		exit(0);
	}
	if (1 != PyDict_Check(pforest)) {
		cerr << " >! Can't find forest type error" << endl;
		exit(0);
	}
	//cout << PyList_Size(pforest) << endl;

	i = _pytype_to_int(PyDict_GetItemString(pforest, "num_x_"),
			" >! error on num_x_");
	j = _pytype_to_int(PyDict_GetItemString(pforest, "num_y_"),
			" >! error on num_y_");
	//i = _pytype_to_int(PyList_GetItem(pforest, 0), " >! error on i");
	//j = _pytype_to_int(PyList_GetItem(pforest, 1), " >! error on j");
	ox = _pytype_to_float(
			PyList_GetItem(PyDict_GetItemString(pforest, "origin_point_"), 0),
			" >! error on i");
	oy = _pytype_to_float(
			PyList_GetItem(PyDict_GetItemString(pforest, "origin_point_"), 1),
			" >! error on j");
	//ox = _pytype_to_float(PyList_GetItem(pforest, 2), " >! error on ox");
	//oy = _pytype_to_float(PyList_GetItem(pforest, 3), " >! error on oy");
	length = _pytype_to_float(PyDict_GetItemString(pforest, "lenght_tree_"),
			" >! error on length tree");
	maxl = _pytype_to_int(PyDict_GetItemString(pforest, "max_level_"),
			" >! error on i");
	if (PyList_Size(pforest) > 6) {
		k = _pytype_to_int(PyDict_GetItemString(pforest, "num_z_"),
				" >! error on k");
		oz = _pytype_to_float(
				PyList_GetItem(PyDict_GetItemString(pforest, "origin_point_"),
						2), " >! error on i");
	}
}

void py_main() {
	PyObject* pParse = NULL;
	Py_Initialize();        //Initialize
	char *path = (char *) "./py/";
	PySys_SetPath(path);
	pParse = PyImport_ImportModule("Larus_config_parse");

	if (!pParse) {
		cerr << " >! Can't find Larus_config_parse " << endl;
		exit(0);
	}
	PyObject* pFunMain = PyObject_GetAttrString(pParse, "main");
	if (!(pFunMain && PyCallable_Check(pFunMain))) {
		cerr << " >! Get function main() in parse is failed!" << endl;
		exit(0);
	} else {
		PyObject* pvreturn = PyObject_CallObject(pFunMain, NULL_PTR);
		if (pvreturn != NULL_PTR) {
			if (0 == PyInt_AsLong(pvreturn)) {
				Py_DECREF(pvreturn);
				exit(0);
			}
		}
	}

	//
	//set tree according to _complete_trees_

	if (get_complete_dim(pParse) == 2) {
		main_2D(pParse);
	}
	cout << "-------so good so far ----------" << endl;
}
template<class FOREST>
void py_set_enable_forest(PyObject* model,  //
		FOREST* f) {
	PyObject* pcts = PyObject_GetAttrString(model, "_complete_trees_");
	if (1 != PyList_Check(pcts)) {
		cerr << " >! Can't find _complete_trees_ : type error" << endl;
		exit(0);
	}
	//cout << "_complete_trees_ " << PyList_Size(pcts) << endl;
	arrayListT<bool> lb(PyList_Size(pcts));
	lb.assign(true);
	for (Py_ssize_t i = 0; i < PyList_Size(pcts); ++i) {
		PyObject* t = PyList_GetItem(pcts, i);
		if (1 != PyDict_Check(t)) {
			cerr << " >! The tree term is not dict : type error" << endl;
			exit(0);
		}
		PyObject* is_enble = PyDict_GetItemString(t, "is_enable");
		if (1 != PyBool_Check(is_enble)) {
			cerr << " >! The tree term is not dict : type error" << endl;
			exit(0);
		}
		lb[int(i)] = (is_enble == Py_True);
	}

	for (arrayListT<bool>::size_type i = 0; i < lb.size(); ++i) {
		f->set_attribution(i, lb[i]);
	}

	//PyEval_CallObject(pFunMain, NULL);
}

Float pFun_call_pyfun(utPointer utp, Float x, Float y, Float z) {
	PyObject* pyfun = CAST(PyObject*, utp);
	PyObject* res = NULL_PTR;
	if (PyFunction_Check(pyfun)) {
		res = PyObject_CallObject(pyfun, Py_BuildValue("ddd", x, y, z));
		return _pytype_to_float(res, " >! return value error ");
	} else {
		cerr << " >! call pyfun uncallable ";
		exit(0);
	}
}

template<class FOREST>
void py_set_variable(
		//
		PyObject* model,  //
		FOREST* f,  //
		arrayListT<string>& lname, arrayList_int& lidx,
		arrayListT<PyObject*>& lpfun) {
	PyObject* pcts = PyObject_GetAttrString(model, "_complete_varlist_");
	if (pcts == NULL_PTR || 1 != PyList_Check(pcts)) {
		cerr << " >! Can't find _complete_variable_ : type error" << endl;
		exit(0);
	}
	int num_v = PyList_Size(pcts);
	lname.reconstruct(num_v);
	lidx.reconstruct(num_v);
	lpfun.reconstruct(num_v);
	for (Py_ssize_t i = 0; i < num_v; ++i) {
		PyObject* t = PyList_GetItem(pcts, i);
		if (1 != PyDict_Check(t)) {
			cerr << " >! The tree term is not dict : type error" << endl;
			exit(0);
		}
		PyObject* name = PyDict_GetItemString(t, "name_");
		if (1 != PyString_Check(name)) {
			cerr << " >! The name out is not string : type error" << endl;
			exit(0);
		}
		lname[i] = _pytype_to_string(name,
				" >! The  is not a string : type error");
		PyObject* idx = PyDict_GetItemString(t, "idx_");
		if (1 != PyInt_Check(idx)) {
			cerr << " >! The tree term is not int : type error" << endl;
			exit(0);
		}
		lidx[i] = _pytype_to_int(idx, " >! The name is not a int : type error");
		PyObject* fun = PyDict_GetItemString(t, "function_");
		if (1 != PyCallable_Check(fun)) {
			cerr << " >! The tree term is not a function : type error" << endl;
			exit(0);
		}
		lpfun[i] = fun;
	}
	//check variable idx
	int max_idx = lidx.findMax();
	for (int i = 0; i < lidx.size(); ++i) {
		if (lidx[i] == 0) {
			cerr << " >! The idx 0 is reserved for index \n";
			exit(0);
		}
	}
	//warning max idx and the lenght of lidx

	//set =========
	resize_array_on_center_leaf(*f, max_idx + 1);
	set_index_on_center_leaf(*f, Idx_IDX);
	//call initial function
	for (int i = 0; i < lidx.size(); ++i) {
		set_value_function_on_leaf( //
				*f, //
				lidx(i), //
				pFun_call_pyfun, //
				lpfun[i]);
	}

	//boundry condition

}

string get_bcdict_key(int idx, int idy, int idz, int dire, int vidx) {
	stringstream ss;
	ss << idx << idy << idz << dire << vidx;
	return ss.str();
}

int fun_set_bc(const BoundaryCondition<Dimension_2D>& bc, pQTNode pn) {
	PyObject* pyfun = CAST(PyObject*, bc.utp);
	PyObject* res = NULL_PTR;
	Float x = pn->cell->getCPX();
	Float y = pn->cell->getCPY();
	Float z = pn->cell->getCPZ();
	if (PyFunction_Check(pyfun)) {
		res = PyObject_CallObject(pyfun, Py_BuildValue("ddd", x, y, z));
		Float resv = _pytype_to_float(res, " >! return value error ");
		pn->data->aCenterData[bc.value_idx] = resv;
	} else {
		cerr << " >! call pyfun uncallable ";
		exit(0);
	}
	return 1;
}

template<class BCM>
void py_set_boundary( //
		PyObject* model,  //
		BCM* bcm, const arrayList_int& lidx) {
	PyObject* pcts = PyObject_GetAttrString(model, "_complete_bcdict_");
	if (pcts == NULL_PTR || 1 != PyDict_Check(pcts)) {
		cerr << " >! Can't find _complete_bcdict_ : type error" << endl;
		exit(0);
	}

	PyObject* pobc;
	PyObject* pofun_item;
	PyObject* pofun;
	for (int i = 0; i < bcm->pforest->iLen(); ++i) {
		for (int j = 0; j < bcm->pforest->jLen(); ++j) {
			for (int k = 0;
					k
							< ((bcm->pforest->get_dim() == 3) ?
									bcm->pforest->kLen() : 1); ++k) {
				pQuadTree pt = bcm->pforest->getpTree(i, j, k);

				if (pt != NULL_PTR) {
					//direction
					for (int ii = 4;
							ii <= ((bcm->pforest->get_dim() == 3) ? 9 : 7);
							ii++) {
						for (int iv = 0; iv < lidx.size(); ++iv) {  //variable
							if (pt->getNeighborpTree(toDirection(ii)) == NULL) {
								//cout<<"in "<< i<<" "<< j<<" "<< k <<"d "<<ii<<" iv "<<lidx[iv]<<endl;
								string key = get_bcdict_key(i, j, k, ii,
										lidx[iv]);
								//cout<<"key "<<key<<endl;
								pobc = PyDict_GetItemString(pcts, key.c_str());
								pofun_item = PyDict_GetItemString(pobc,
										"bc_fun");
								pofun = PyTuple_GetItem(pofun_item, 2);
								BoundaryCondition<Dimension_2D> bc;
								bc.direction = toDirection(ii);
								bc.tree_idx = bcm->pforest->to_1d_idx(i, j, k);
								bc.value_idx = lidx[iv];
								bc.utp = pofun;
								bc.pfun = fun_set_bc;
								bcm->add_BC(bc);
							}
						}
					}
				}
			}
		}
	}
	bcm->set_bc();
	//bcm->show_set_of_boundary_condition();
	cout << "end of py set boundary ------------\n";
}

int _find_idx(const string name, const arrayListT<string>& lname,
		const arrayList_int& lidx) {
	int idx = 0;
	for (int i = 0; i < lname.size(); i++) {
		if (lname[i] == name) {
			idx = i;
			return lidx[idx];
		}
	}
	assert(false);
	return 0;
}

Float pfunres(Float x, Float y, Float z) {
	return y * (1 - y) * x * x * x;
}

void main_2D(PyObject* pParse) {
	cout << ">  process to main2D" << endl;
	int ni, nj, nk, maxl;
	Float ox, oy, oz, length;
	py_construct_param_forest(pParse, ni, nj, ox, oy, length, maxl, nk, oz);
	Forest2D forest(ni, nj, ox, oy, length, maxl);
	py_set_enable_forest(pParse, &forest);
	for (int i = 0; i < forest.size(); i++) {
		pQuadTree tree = forest.getpTree_1d(i);
		if (tree != NULL_PTR) {
			tree->CreatFullTree();
		}
	}
	forest.ConnectTrees();
//read and set variable _complete_varlist_
	arrayListT<string> lname;
	arrayList_int lidx;
	arrayListT<PyObject*> lpfun;
	py_set_variable(pParse, &forest, lname, lidx, lpfun);
	forest.show_info();
	forest.show();
	BCManager<Dimension_2D> bcm(&forest);
	bcm.new_ghost_nodes();
	py_set_boundary(pParse, &bcm, lidx);
	//bcm.show_tree_boundary();

	int idx_beta = _find_idx("beta", lname, lidx);
	int idx_phi = _find_idx("phi", lname, lidx);
	int idx_f = _find_idx("f", lname, lidx);
	cout << "bata idx " << idx_beta << "\n";
	cout << "phi  idx " << idx_phi << "\n";
	cout << "f    idx " << idx_f << "\n";
	Poisson_Eq<Dimension_2D> pe(&forest, &bcm, idx_beta, idx_phi, idx_f);
	bcm.show_boundary_contour(idx_phi);
	Float tol = 1e-6;
	//pe.set_f_term(f_fun);
	pe.solve(tol);
	pe.show();
	//the forest for compare ---------------------
	//Forest2D forest2(ni, nj, ox, oy, length, maxl);
	//for (int i = 0; i < forest2.size(); i++) {
	//	pQuadTree tree = forest2.getpTree_1d(i);
	//	if (tree != NULL_PTR) {
	//		tree->CreatFullTree();
	//	}
	//}
	//forest2.ConnectTrees();
	//forest2.show_info();
	//resize_array_on_center_leaf(forest2, 2);
	//set_index_on_center_leaf(forest2, Idx_IDX);
	//set_scalar_on_leaf_by_function(forest2, 1, pfunres);
	//show_gnuplot_as_contour(forest2, 1);
}

}

#endif /* TEST_PYTHON_H_ */
