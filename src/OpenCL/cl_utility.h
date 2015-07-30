/*
 * cl_utility.h
 *
 *  Created on: Jul 3, 2015
 *      Author: zhou
 */

#ifndef CL_UTILITY_H_
#define CL_UTILITY_H_

#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#ifdef __cplusplus
#include <cstdio>
#endif

namespace cl {

const char *err_code(cl_int err_in);

void check_error(cl_int err, const std::string& operation, const std::string& filename, int line);

#define checkError(E, S) check_error(E,S,__FILE__,__LINE__)

int show_opencl_info();

char* getKernelSource(const std::string& filename);

}

#endif /* OPENCL_CL_UTILITY_H_ */
