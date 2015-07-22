/************************
 //  \file   IO_tablefile.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   2 sept. 2014 
 ***********************/

#include "IO_tablefile.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>

using namespace std;
namespace Larus {
string TableFile::getFilename() {
	return _filename;
}

int TableFile::getRow() {
	return val.iLen();
}
int TableFile::getCol() {
	return val.jLen();
}
int TableFile::countIdx() {
	return idx.size();
}
//Float getIdx(int i);
//string getIdxName(int i);

TableFile::TableFile() {
	_filename = "";
	idx = arrayList();
	idx_name = arrayListT<string>();
	row_name = arrayListT<string>();
	col_name = arrayListT<string>();
	val = Matrix();
}
TableFile::TableFile(string name) {
	_filename = name;
	idx = arrayList();
	idx_name = arrayListT<string>();
	row_name = arrayListT<string>();
	col_name = arrayListT<string>();
	val = Matrix();
}
TableFile::TableFile(string name, Matrix &m) {
	_filename = name;
	idx = arrayList();
	idx_name = arrayListT<string>();
	row_name = arrayListT<string>(m.jLen());
	col_name = arrayListT<string>(m.jLen());
	val = m;
}

void TableFile::setRowName(idx_t i, string name) {
	if (i < 0 || i >= row_name.size()) {
		cerr << ">! set Row name error, i= " << i << "out of range";
	} else {
		row_name[i] = name;
	}
}
void TableFile::setColName(idx_t j, string name) {
	if (j < 0 || j >= col_name.size()) {
		cerr << ">! set Col name error, j= " << j << "out of range";
	} else {
		col_name[j] = name;
	}
}

void TableFile::read() {
	ifstream ifs;
	ifs.open(this->_filename.c_str(), std::ifstream::in);
	if (!ifs.is_open()) {
		cerr << "Error opening file \n";
		exit(0);
	} else {
		idx_t iline=0;
		idx_t r=0,c=0;
		while (!ifs.eof()) {
			string sline;
			getline(ifs, sline, '\n');
			if (sline.length() == 0) {
				//cout<<"Empty"<<endl;
				//emptyline++;
			} else {
				istringstream istr(sline);
				string s;
				istr >> s;
				if (s == "#") {
					istr>>s;
					if(s=="NAME:"){

					}else if(s=="ROW:"){
						istr>>r;
					}else if(s=="COL:"){
						istr>>c;
					}
				}else{
					istr.seekg(0,istr.beg);
					if(iline==0){
						this->val.reconstruct(r,c);
					}
					for(int j=0;j<this->val.jLen();j++){
						istr >> this->val[iline][j];
					}
					iline++;
				}
				//mgroup.appendRow(a);
				//cout <<a[0]<<"# #"<<a[1]<<"# #"<<a[2]<<endl;
				//segline++;
			}
		}
	}
}

void TableFile::output() const {
	FILE *data;
	data = fopen(_filename.c_str(), "w"); //write

	fprintf(data, "# NAME: %s\n", _filename.c_str());
	fprintf(data, "# ROW: %d\n", val.iLen());
	fprintf(data, "# COL: %d\n", val.jLen());
	for (int i = 0; i < val.iLen(); i++) {
		for (int j = 0; j < val.jLen(); j++) {
			fprintf(data, "%lf ", val[i][j]);
		}
		fprintf(data, "\n");
	}
	fclose(data);
}

}
