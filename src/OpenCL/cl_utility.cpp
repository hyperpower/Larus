#include "cl_utility.h"
#include <string>

namespace cl {

const char *err_code(cl_int err_in) {
	switch (err_in) {
	case CL_SUCCESS:
		return (char*) "CL_SUCCESS";
	case CL_DEVICE_NOT_FOUND:
		return (char*) "CL_DEVICE_NOT_FOUND";
	case CL_DEVICE_NOT_AVAILABLE:
		return (char*) "CL_DEVICE_NOT_AVAILABLE";
	case CL_COMPILER_NOT_AVAILABLE:
		return (char*) "CL_COMPILER_NOT_AVAILABLE";
	case CL_MEM_OBJECT_ALLOCATION_FAILURE:
		return (char*) "CL_MEM_OBJECT_ALLOCATION_FAILURE";
	case CL_OUT_OF_RESOURCES:
		return (char*) "CL_OUT_OF_RESOURCES";
	case CL_OUT_OF_HOST_MEMORY:
		return (char*) "CL_OUT_OF_HOST_MEMORY";
	case CL_PROFILING_INFO_NOT_AVAILABLE:
		return (char*) "CL_PROFILING_INFO_NOT_AVAILABLE";
	case CL_MEM_COPY_OVERLAP:
		return (char*) "CL_MEM_COPY_OVERLAP";
	case CL_IMAGE_FORMAT_MISMATCH:
		return (char*) "CL_IMAGE_FORMAT_MISMATCH";
	case CL_IMAGE_FORMAT_NOT_SUPPORTED:
		return (char*) "CL_IMAGE_FORMAT_NOT_SUPPORTED";
	case CL_BUILD_PROGRAM_FAILURE:
		return (char*) "CL_BUILD_PROGRAM_FAILURE";
	case CL_MAP_FAILURE:
		return (char*) "CL_MAP_FAILURE";
		//case CL_MISALIGNED_SUB_BUFFER_OFFSET:
		//    return (char*)"CL_MISALIGNED_SUB_BUFFER_OFFSET";
		//case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
		//    return (char*)"CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
	case CL_INVALID_VALUE:
		return (char*) "CL_INVALID_VALUE";
	case CL_INVALID_DEVICE_TYPE:
		return (char*) "CL_INVALID_DEVICE_TYPE";
	case CL_INVALID_PLATFORM:
		return (char*) "CL_INVALID_PLATFORM";
	case CL_INVALID_DEVICE:
		return (char*) "CL_INVALID_DEVICE";
	case CL_INVALID_CONTEXT:
		return (char*) "CL_INVALID_CONTEXT";
	case CL_INVALID_QUEUE_PROPERTIES:
		return (char*) "CL_INVALID_QUEUE_PROPERTIES";
	case CL_INVALID_COMMAND_QUEUE:
		return (char*) "CL_INVALID_COMMAND_QUEUE";
	case CL_INVALID_HOST_PTR:
		return (char*) "CL_INVALID_HOST_PTR";
	case CL_INVALID_MEM_OBJECT:
		return (char*) "CL_INVALID_MEM_OBJECT";
	case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
		return (char*) "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
	case CL_INVALID_IMAGE_SIZE:
		return (char*) "CL_INVALID_IMAGE_SIZE";
	case CL_INVALID_SAMPLER:
		return (char*) "CL_INVALID_SAMPLER";
	case CL_INVALID_BINARY:
		return (char*) "CL_INVALID_BINARY";
	case CL_INVALID_BUILD_OPTIONS:
		return (char*) "CL_INVALID_BUILD_OPTIONS";
	case CL_INVALID_PROGRAM:
		return (char*) "CL_INVALID_PROGRAM";
	case CL_INVALID_PROGRAM_EXECUTABLE:
		return (char*) "CL_INVALID_PROGRAM_EXECUTABLE";
	case CL_INVALID_KERNEL_NAME:
		return (char*) "CL_INVALID_KERNEL_NAME";
	case CL_INVALID_KERNEL_DEFINITION:
		return (char*) "CL_INVALID_KERNEL_DEFINITION";
	case CL_INVALID_KERNEL:
		return (char*) "CL_INVALID_KERNEL";
	case CL_INVALID_ARG_INDEX:
		return (char*) "CL_INVALID_ARG_INDEX";
	case CL_INVALID_ARG_VALUE:
		return (char*) "CL_INVALID_ARG_VALUE";
	case CL_INVALID_ARG_SIZE:
		return (char*) "CL_INVALID_ARG_SIZE";
	case CL_INVALID_KERNEL_ARGS:
		return (char*) "CL_INVALID_KERNEL_ARGS";
	case CL_INVALID_WORK_DIMENSION:
		return (char*) "CL_INVALID_WORK_DIMENSION";
	case CL_INVALID_WORK_GROUP_SIZE:
		return (char*) "CL_INVALID_WORK_GROUP_SIZE";
	case CL_INVALID_WORK_ITEM_SIZE:
		return (char*) "CL_INVALID_WORK_ITEM_SIZE";
	case CL_INVALID_GLOBAL_OFFSET:
		return (char*) "CL_INVALID_GLOBAL_OFFSET";
	case CL_INVALID_EVENT_WAIT_LIST:
		return (char*) "CL_INVALID_EVENT_WAIT_LIST";
	case CL_INVALID_EVENT:
		return (char*) "CL_INVALID_EVENT";
	case CL_INVALID_OPERATION:
		return (char*) "CL_INVALID_OPERATION";
	case CL_INVALID_GL_OBJECT:
		return (char*) "CL_INVALID_GL_OBJECT";
	case CL_INVALID_BUFFER_SIZE:
		return (char*) "CL_INVALID_BUFFER_SIZE";
	case CL_INVALID_MIP_LEVEL:
		return (char*) "CL_INVALID_MIP_LEVEL";
	case CL_INVALID_GLOBAL_WORK_SIZE:
		return (char*) "CL_INVALID_GLOBAL_WORK_SIZE";
		//case CL_INVALID_PROPERTY:
		//    return (char*)"CL_INVALID_PROPERTY";

	default:
		return (char*) "UNKNOWN ERROR";
	}
}



void check_error(cl_int err, const std::string& operation, const std::string& filename, int line)
{

    if (err != CL_SUCCESS)
    {
        fprintf(stderr, "Error during operation '%s', ", operation.c_str());
        fprintf(stderr, "in '%s' on line %d\n", filename.c_str(), line);
        fprintf(stderr, "Error code was \"%s\" (%d)\n", err_code(err), err);
        exit(EXIT_FAILURE);
    }
}


int show_opencl_info() {
	cl_int err;
	// Find the number of OpenCL platforms
	cl_uint num_platforms;
	err = clGetPlatformIDs(0, NULL, &num_platforms);
	checkError(err, "Finding platforms");
	if (num_platforms == 0) {
		printf("Found 0 platforms!\n");
		return EXIT_FAILURE;
	}
	// Create a list of platform IDs
	cl_platform_id platform[num_platforms];
	err = clGetPlatformIDs(num_platforms, platform, NULL);
	checkError(err, "Getting platforms");

	printf("\nNumber of OpenCL platforms: %d\n", num_platforms);
	printf("\n-------------------------\n");

	// Investigate each platform
	for (cl_uint i = 0; i < num_platforms; i++) {
		cl_char string[10240] = { 0 };
		// Print out the platform name
		err = clGetPlatformInfo(platform[i], CL_PLATFORM_NAME, sizeof(string),
				&string, NULL);
		checkError(err, "Getting platform name");
		printf("Platform: %s\n", string);

		// Print out the platform vendor
		err = clGetPlatformInfo(platform[i], CL_PLATFORM_VENDOR, sizeof(string),
				&string, NULL);
		checkError(err, "Getting platform vendor");
		printf("Vendor: %s\n", string);

		// Print out the platform OpenCL version
		err = clGetPlatformInfo(platform[i], CL_PLATFORM_VERSION,
				sizeof(string), &string, NULL);
		checkError(err, "Getting platform OpenCL version");
		printf("Version: %s\n", string);

		// Count the number of devices in the platform
		cl_uint num_devices;
		err = clGetDeviceIDs(platform[i], CL_DEVICE_TYPE_ALL, 0, NULL,
				&num_devices);
		checkError(err, "Finding devices");

		// Get the device IDs
		cl_device_id device[num_devices];
		err = clGetDeviceIDs(platform[i], CL_DEVICE_TYPE_ALL, num_devices,
				device, NULL);
		checkError(err, "Getting devices");
		printf("Number of devices: %d\n", num_devices);

		// Investigate each device
		for (cl_uint j = 0; j < num_devices; j++) {
			printf("\t-------------------------\n");

			// Get device name
			err = clGetDeviceInfo(device[j], CL_DEVICE_NAME, sizeof(string),
					&string, NULL);
			checkError(err, "Getting device name");
			printf("\t Name: %s\n", string);

			// Get device OpenCL version
			//err = clGetDeviceInfo(device[j], CL_DEVICE_VERSION , sizeof(string), &string, NULL);
			checkError(err, "Getting device OpenCL C version");
			printf("\t Version: %s\n", string);

			// Get Max. Compute units
			cl_uint num;
			err = clGetDeviceInfo(device[j], CL_DEVICE_MAX_COMPUTE_UNITS,
					sizeof(cl_uint), &num, NULL);
			checkError(err, "Getting device max compute units");
			printf("\t Max. Compute Units: %d\n", num);

			// Get local memory size
			cl_ulong mem_size;
			err = clGetDeviceInfo(device[j], CL_DEVICE_LOCAL_MEM_SIZE,
					sizeof(cl_ulong), &mem_size, NULL);
			checkError(err, "Getting device local memory size");
			printf("\t Local Memory Size: %lu KB\n", mem_size / 1024);

			// Get global memory size
			err = clGetDeviceInfo(device[j], CL_DEVICE_GLOBAL_MEM_SIZE,
					sizeof(cl_ulong), &mem_size, NULL);
			checkError(err, "Getting device global memory size");
			printf("\t Global Memory Size: %lu MB\n",
					mem_size / (1024 * 1024));

			// Get maximum buffer alloc. size
			err = clGetDeviceInfo(device[j], CL_DEVICE_MAX_MEM_ALLOC_SIZE,
					sizeof(cl_ulong), &mem_size, NULL);
			checkError(err, "Getting device max allocation size");
			printf("\t Max Alloc Size: %lu MB\n", mem_size / (1024 * 1024));
			int num_double = mem_size/sizeof(double);
			int num_float = mem_size/sizeof(float);
			printf("\t Max Alloc Num double: %d\n", num_double);
			printf("\t Max Alloc Num float: %d\n", num_float);

			// Get work-group size information
			size_t size;
			err = clGetDeviceInfo(device[j], CL_DEVICE_MAX_WORK_GROUP_SIZE,
					sizeof(size_t), &size, NULL);
			checkError(err, "Getting device max work-group size");
			printf("\t Max Work-group Total Size: %ld\n", size);

			// Find the maximum dimensions of the work-groups
			err = clGetDeviceInfo(device[j], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,
					sizeof(cl_uint), &num, NULL);
			checkError(err, "Getting device max work-item dims");
			// Get the max. dimensions of the work-groups
			size_t dims[num];
			err = clGetDeviceInfo(device[j], CL_DEVICE_MAX_WORK_ITEM_SIZES,
					sizeof(dims), &dims, NULL);
			checkError(err, "Getting device max work-item sizes");
			printf("\t Max Work-group Dims: ( ");
			for (size_t k = 0; k < num; k++) {
				printf("%ld ", dims[k]);
			}
			printf(")\n");

			printf("\t-------------------------\n");
		}

		printf("\n-------------------------\n");
	}

	return EXIT_SUCCESS;
}

char * getKernelSource(const std::string& filename)
{

    FILE *file = fopen(filename.c_str(), "r");
    if (!file)
    {
        fprintf(stderr, "Error: Could not open kernel source file\n");
        exit(EXIT_FAILURE);
    }
    fseek(file, 0, SEEK_END);
    int len = ftell(file) + 1;
    rewind(file);

    char *source = (char *)calloc(sizeof(char), len);
    if (!source)
    {
        fprintf(stderr, "Error: Could not allocate memory for source string\n");
        exit(EXIT_FAILURE);
    }
    fread(source, sizeof(char), len, file);
    fclose(file);
    return source;
}

}
