/************************
 //  \file   IO.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   3 déc. 2014 
 ***********************/

#include "IO.h"
#include <string>

using namespace std;

namespace Larus {
FILE* open_file(string filename, int mode) {
	FILE *data = NULL;
	if (mode == 1) {
		data = fopen(filename.c_str(), "w"); //write
	} else if (mode == 2) {
		data = fopen(filename.c_str(), "a"); //append
		if (data == NULL) {
			std::cerr << "!> Can't find file. " << filename.c_str() << " \n";
			exit(-1);
		}
	}
	if (data == NULL) {
		std::cerr << "!> Open file error! " << filename.c_str() << " \n";
		exit(-1);
	}
	return data;
}

void creat_empty_file(string filename) {
	FILE *data;
	data = fopen(filename.c_str(), "w"); //write
	fclose(data);
}

void append_empty_line(string filename){
	FILE* data = open_file(filename, 2); //append
	fprintf(data,"\n");
	fclose(data);
}

void open_file(ofstream& fs, string filename, int mode) {
	if (mode == 1) {
		fs.open (filename.c_str(), fstream::out ); //write
	} else if (mode == 2) {
		fs.open (filename.c_str(), fstream::out | fstream::app); //append
		if (!fs.is_open()) {
			std::cerr << "!> Can't find file. " << filename.c_str() << " \n";
			exit(-1);
		}
	}
	if (!fs.is_open()) {
		std::cerr << "!> Open file error! " << filename.c_str() << " \n";
		exit(-1);
	}
}

bool file_access_check(const std::string &filename, int mode)
{
	ASSERT(mode >= 0 && mode <= 7);
	//  int _access(const char *path, int mode);
	//  returns 0 if the file has the given mode,
	//  it returns -1 if the named file does not exist or is not accessible in
	//  the given mode
	// mode = 0 (F_OK) (default): checks file for existence only
	// mode = 1 (X_OK): execution permission
	// mode = 2 (W_OK): write permission
	// mode = 4 (R_OK): read permission
	// mode = 6       : read and write permission
	// mode = 7       : read, write and execution permission
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
	if (_access(filename.c_str(), mode) == 0)
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
	if (access(filename.c_str(), mode) == 0)
#endif
			{
		return true;
	} else {
		return false;
	}
}

}
