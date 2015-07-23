/*
 * cl_utility.h
 *
 *  Created on: Jul 3, 2015
 *      Author: zhou
 */

#ifndef CL_UTILITY_H_
#define CL_UTILITY_H_

#include "cl.hpp"

namespace cl {

#include <cstdio>
#include <cstdlib>
#include <iostream>

const char * helloStr = "__kernel void "
		"hello(void) "
		"{ "
		"  "
		"} ";

inline cl_int show_device_info() {
	cl_int err_code = CL_SUCCESS;
	try {
		std::vector<Platform> platforms;
		Platform::get(&platforms);
		cout << "Platform size :" << platforms.size() << "\n";
		if (platforms.size() == 0) {
			cout << "Platform size 0\n";
			return -1;
		}
		// Investigate each platform
		int pf_idx = 0;
		for (std::vector<Platform>::iterator iter = platforms.begin();
				iter != platforms.end(); ++iter) {
			std::cout << "Platforms id:" << pf_idx << "\n";
			// Print out the platform name
			std::string name;
			iter->getInfo(CL_PLATFORM_NAME, &name);
			std::cout << " > Name       " << name << "\n";
			// Print out the platform vendor
			std::string vendor;
			iter->getInfo(CL_PLATFORM_VENDOR, &vendor);
			std::cout << " > Vendor     " << name << "\n";
			// Print out the platform version
			std::string version;
			iter->getInfo(CL_PLATFORM_VERSION, &version);
			std::cout << " > Version    " << version << "\n";
			//
			std::vector<Device> devices;
			iter->getDevices(CL_DEVICE_TYPE_ALL, &devices);
			std::cout << "+> Num Devices " << devices.size() << "\n";

			int de_idx = 0;
			for (std::vector<Device>::iterator iter = devices.begin();
					iter != devices.end(); ++iter) {
				std::cout << "  +> " << de_idx << "--------\n";
				std::string str;
				iter->getInfo(CL_DEVICE_NAME, &str);
				std::cout << "   > Name               : " << str << "\n";
				iter->getInfo(CL_DEVICE_VERSION, &str);
				std::cout << "   > Version            : " << str << "\n";
				cl_uint num;
				iter->getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &num);
				std::cout << "   > Max Compute Unit   : " << num << "\n";
				cl_ulong mem_size;
				iter->getInfo(CL_DEVICE_LOCAL_MEM_SIZE, &mem_size);
				std::cout << "   > Local Memory Size  : " << mem_size / 1024
						<< " KB\n";
				iter->getInfo(CL_DEVICE_GLOBAL_MEM_SIZE, &mem_size);
				std::cout << "   > Global Memory Size : "
						<< mem_size / 1024 / 1024 << " MB\n";
				iter->getInfo(CL_DEVICE_MAX_MEM_ALLOC_SIZE, &mem_size);
				std::cout << "   > Max Allocation Size: "
						<< mem_size / 1024 / 1024 << " MB\n";
				std::size_t s;
				iter->getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &s);
				std::cout << "   > Max Work-group Size: " << s << "\n";

				iter->getInfo(CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, &num);
				std::cout << "   > Max Work-group Dim : " << num << "\n";
				std::vector<std::size_t> dims;
				iter->getInfo(CL_DEVICE_MAX_WORK_ITEM_SIZES, &dims);
				std::cout << "   > Work-group Dim Size: ";
				for (std::size_t k = 0; k < num; k++) {
					std::cout << dims[k] << " ";
				}
				std::cout << "\n";
				de_idx++;
			}
			pf_idx++;
		}
	} catch (Error& err) {
		std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")"
				<< std::endl;
	}
	return err_code;
}

inline cl_int work_flow_control() {
	cl_int err_code = CL_SUCCESS;
	try {
		cl_uint pf_Idx = 0;
		cl_uint de_Idx = 1;

		std::vector<Platform> platforms;
		Platform::get(&platforms);
		if (platforms.size() == 0) {
			return CL_DEVICE_NOT_FOUND;
		}
		// get devices
		std::vector<Device> devices;
		platforms[pf_Idx].getDevices(CL_DEVICE_TYPE_ALL, &devices);
		//
		if (de_Idx >= devices.size()) {
			std::cerr << " de_Idx out of range\n";
			return CL_DEVICE_NOT_FOUND;
		}
		Device& de = devices[de_Idx];
		// Create a compute context
		Context context(devices, NULL, NULL, &err_code);
		// Create a command queue
		CommandQueue cq(context, de, 0, &err_code);

	} catch (Error& err) {
		std::cerr << "ERROR: " << err.what() << "(" << err.err() << ")"
				<< std::endl;
	}
	return err_code;
}

}

#endif /* OPENCL_CL_UTILITY_H_ */
