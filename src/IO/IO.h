/************************
 //  \file   IO.h
 //  \brief
 // 
 //  \author czhou
 //  \date   3 déc. 2014 
 ***********************/
#ifndef IO_H_
#define IO_H_

#include <iostream>
#include <fstream>
#include <unistd.h>

#include "../TypeDef.h"
#include <string>
#include <fstream>


namespace Larus{
	//c-style
    void creat_empty_file(std::string filename);
    void append_empty_line(std::string filename);
	FILE* open_file(std::string, int);
	//c++style
	void open_file(std::ofstream& , std::string, int);
	bool file_access_check(const std::string &filename, int=1);
}



#endif /* IO_H_ */
