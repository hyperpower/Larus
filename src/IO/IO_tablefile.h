/************************
 //  \file   IO_table.h
 //  \brief
 // 
 //  \author zhou
 //  \date   22 août 2014 
 ***********************/
#ifndef IO_TABLEFILE_H_
#define IO_TABLEFILE_H_

#include "../TypeDef.h"
#include "../Utility/Array.h"
#include "../Algebra/Matrix.h"
#include <string>

using namespace std;
namespace Larus{

class TableFile{
private:
	string _filename;
public:
	arrayList idx;
	arrayListT<string> idx_name;
	arrayListT<string> row_name;
	arrayListT<string> col_name;
	Matrix val;

	//Function
	string getFilename();
	int getRow();
	int getCol();
	int countIdx();
	//Float getIdx(int i);
	//string getIdxName(int i);
	void read();
	TableFile();
	TableFile(string name);
	TableFile(string name, Matrix &m);
	void setRowName(idx_t i, string name);
	void setColName(idx_t j, string name);

	void output() const;

};


}




#endif /* IO_TABLE_H_ */
