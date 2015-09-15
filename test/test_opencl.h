#ifndef _TEST_OPENCL_H_
#define _TEST_OPENCL_H_

#include "../src/Algebra/Arithmetic.h"
#include "../src/Utility/ArrayList.h"
#include "../src/OpenCL/cl_utility.h"

using namespace Larus;

namespace cl {

#define LENGTH (1024*1024)
#define TOL    (0.001)
typedef double vt;

int vector_add() {

	int err;               // error code returned from OpenCL calls

	vt* h_a = (vt*) calloc(LENGTH, sizeof(vt));       // a vector
	vt* h_b = (vt*) calloc(LENGTH, sizeof(vt));       // b vector
	vt* h_c = (vt*) calloc(LENGTH, sizeof(vt)); // c vector (a+b) returned from the compute device

	unsigned int correct;           // number of correct results

	size_t global;                  // global domain size

	cl_context context;       // compute context
	cl_command_queue commands;      // compute command queue
	cl_program program;       // compute program
	cl_kernel ko_vadd;       // compute kernel

	cl_mem d_a;                    // device memory used for the input  a vector
	cl_mem d_b;                    // device memory used for the input  b vector
	cl_mem d_c;                    // device memory used for the output c vector

	// Fill vectors a and b with random float values
	int i = 0;
	int count = LENGTH;
	for (i = 0; i < count; i++) {
		h_a[i] = rand() / (vt) RAND_MAX;
		h_b[i] = rand() / (vt) RAND_MAX;
	}

	// Set up platform and GPU device

	cl_uint numPlatforms;

	// Find number of platforms
	err = clGetPlatformIDs(0, NULL, &numPlatforms);
	checkError(err, "Finding platforms");
	if (numPlatforms == 0) {
		printf("Found 0 platforms!\n");
		return EXIT_FAILURE;
	}

	// Get all platforms
	cl_platform_id Platform[numPlatforms];
	err = clGetPlatformIDs(numPlatforms, Platform, NULL);
	checkError(err, "Getting platforms");

	cl_uint num_devices;
	err = clGetDeviceIDs(Platform[0], CL_DEVICE_TYPE_ALL, 0, NULL,
			&num_devices);
	checkError(err, "Finding devices");

	// Get the device IDs
	cl_device_id device[num_devices];
	err = clGetDeviceIDs(Platform[0], CL_DEVICE_TYPE_ALL, num_devices, device,
			NULL);
	checkError(err, "Getting devices");

	cl_device_id device_id = device[1];

	// Create a compute context
	context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
	checkError(err, "Creating context");

	// Create a command queue
	commands = clCreateCommandQueue(context, device_id, 0, &err);
	checkError(err, "Creating command queue ");

	char* KernelSource = getKernelSource("src/OpenCL/vector.cl");
	// Create the compute program from the source buffer
	program = clCreateProgramWithSource(context, 1,
			(const char **) &KernelSource, NULL, &err);
	checkError(err, "Creating program ");

	// Build the program
	err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	if (err != CL_SUCCESS) {
		size_t len;
		char buffer[2048];

		printf("Error: Failed to build program executable!\n%s\n",
				err_code(err));
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG,
				sizeof(buffer), buffer, &len);
		printf("%s\n", buffer);
		return EXIT_FAILURE;
	}

	// Create the compute kernel from the program
	ko_vadd = clCreateKernel(program, "vec_add", &err);
	checkError(err, "Creating kernel");

	// Create the input (a, b) and output (c) arrays in device memory
	d_a = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(vt) * count, NULL,
			&err);
	checkError(err, "Creating buffer d_a");

	d_b = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(vt) * count, NULL,
			&err);
	checkError(err, "Creating buffer d_b");

	d_c = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(vt) * count, NULL,
			&err);
	checkError(err, "Creating buffer d_c");

	// Write a and b vectors into compute device memory
	err = clEnqueueWriteBuffer(commands, d_a, CL_TRUE, 0, sizeof(vt) * count,
			h_a, 0, NULL, NULL);
	checkError(err, "Copying h_a to device at d_a");

	err = clEnqueueWriteBuffer(commands, d_b, CL_TRUE, 0, sizeof(vt) * count,
			h_b, 0, NULL, NULL);
	checkError(err, "Copying h_b to device at d_b");

	// Set the arguments to our compute kernel
	err = clSetKernelArg(ko_vadd, 0, sizeof(cl_mem), &d_a);
	err |= clSetKernelArg(ko_vadd, 1, sizeof(cl_mem), &d_b);
	err |= clSetKernelArg(ko_vadd, 2, sizeof(cl_mem), &d_c);
	err |= clSetKernelArg(ko_vadd, 3, sizeof(unsigned int), &count);
	checkError(err, "Setting kernel arguments");

	//double rtime = wtime();

	// Execute the kernel over the entire range of our 1d input data set
	// letting the OpenCL runtime choose the work-group size
	global = count;
	err = clEnqueueNDRangeKernel(commands, ko_vadd, 1, NULL, &global, NULL, 0,
			NULL, NULL);
	checkError(err, "Enqueueing kernel");

	// Wait for the commands to complete before stopping the timer
	err = clFinish(commands);
	checkError(err, "Waiting for kernel to finish");

	//rtime = wtime() - rtime;
	//printf("\nThe kernel ran in %lf seconds\n",rtime);

	// Read back the results from the compute device
	err = clEnqueueReadBuffer(commands, d_c, CL_TRUE, 0, sizeof(vt) * count,
			h_c, 0, NULL, NULL);
	if (err != CL_SUCCESS) {
		printf("Error: Failed to read output array!\n%s\n", err_code(err));
		exit(1);
	}

	// Test the results
	correct = 0;
	vt tmp;

	for (i = 0; i < count; i++) {
		tmp = h_a[i] + h_b[i];     // assign element i of a+b to tmp
		//tmp -= h_c[i];             // compute deviation of expected and output result
		if (tmp * tmp < TOL * TOL) { // correct if square deviation is less than tolerance squared
			correct++;
			printf(" tmp %f h_a %f h_b %f h_c %f \n", tmp, h_a[i], h_b[i],
					h_c[i]);
		} else {
			printf(" tmp %f h_a %f h_b %f h_c %f \n", tmp, h_a[i], h_b[i],
					h_c[i]);
		}
	}

	// summarise results
	printf("C = A+B:  %d out of %d results were correct.\n", correct, count);

	// cleanup then shutdown
	clReleaseMemObject(d_a);
	clReleaseMemObject(d_b);
	clReleaseMemObject(d_c);
	clReleaseProgram(program);
	clReleaseKernel(ko_vadd);
	clReleaseCommandQueue(commands);
	clReleaseContext(context);

	free(h_a);
	free(h_b);
	free(h_c);

	return 0;
}

int vector_dot() {

	int err;               // error code returned from OpenCL calls

	vt* h_a = (vt*) calloc(LENGTH, sizeof(vt));       // a vector
	vt* h_b = (vt*) calloc(LENGTH, sizeof(vt));       // b vector
	vt* h_c = (vt*) calloc(LENGTH, sizeof(vt)); // c vector (a+b) returned from the compute device

	unsigned int correct;           // number of correct results

	size_t global;                  // global domain size

	cl_context context;       // compute context
	cl_command_queue commands;      // compute command queue
	cl_program program;       // compute program
	cl_kernel ko_vecdot[2];  // compute kernel

	cl_mem d_a;                    // device memory used for the input  a vector
	cl_mem d_b;                    // device memory used for the input  b vector
	cl_mem d_c;                    // device memory used for the output c vector

	// Fill vectors a and b with random float values
	int i = 0;
	int count = LENGTH;
	for (i = 0; i < count; i++) {
		h_a[i] = rand() / (vt) RAND_MAX;
		h_b[i] = rand() / (vt) RAND_MAX;
		//h_a[i] = 2;
		//h_b[i] = 2;
	}

	// Set up platform and GPU device

	cl_uint numPlatforms;

	// Find number of platforms
	err = clGetPlatformIDs(0, NULL, &numPlatforms);
	checkError(err, "Finding platforms");
	if (numPlatforms == 0) {
		printf("Found 0 platforms!\n");
		return EXIT_FAILURE;
	}

	// Get all platforms
	cl_platform_id Platform[numPlatforms];
	err = clGetPlatformIDs(numPlatforms, Platform, NULL);
	checkError(err, "Getting platforms");

	cl_uint num_devices;
	err = clGetDeviceIDs(Platform[0], CL_DEVICE_TYPE_ALL, 0, NULL,
			&num_devices);
	checkError(err, "Finding devices");

	// Get the device IDs
	cl_device_id device[num_devices];
	err = clGetDeviceIDs(Platform[0], CL_DEVICE_TYPE_ALL, num_devices, device,
			NULL);
	checkError(err, "Getting devices");

	cl_device_id device_id = device[0];

	// Create a compute context
	context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
	checkError(err, "Creating context");

	// Create a command queue
	commands = clCreateCommandQueue(context, device_id, 0, &err);
	checkError(err, "Creating command queue ");

	char* KernelSource = getKernelSource("src/OpenCL/vector.cl");
	// Create the compute program from the source buffer
	program = clCreateProgramWithSource(context, 1,
			(const char **) &KernelSource, NULL, &err);
	checkError(err, "Creating program ");

	// Build the program
	err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	if (err != CL_SUCCESS) {
		size_t len;
		char buffer[2048];

		printf("Error: Failed to build program executable!\n%s\n",
				err_code(err));
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG,
				sizeof(buffer), buffer, &len);
		printf("%s\n", buffer);
		return EXIT_FAILURE;
	}

	// Create the compute kernel from the program
	ko_vecdot[0] = clCreateKernel(program, "vec_dot", &err);
	checkError(err, "Creating kernel");
	ko_vecdot[1] = clCreateKernel(program, "vec_sum", &err);
	checkError(err, "Creating kernel");

	// Create the input (a, b) and output (c) arrays in device memory
	d_a = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(vt) * count, NULL,
			&err);
	checkError(err, "Creating buffer d_a");

	d_b = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(vt) * count, NULL,
			&err);
	checkError(err, "Creating buffer d_b");

	d_c = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(vt) * count, NULL,
			&err);
	checkError(err, "Creating buffer d_c");

	cl_mem d_res;
	vt h_res = 0;
	d_res = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(vt), NULL,
				&err);
	checkError(err, "Creating buffer d_c");

	// Write a and b vectors into compute device memory
	err = clEnqueueWriteBuffer(commands, d_a, CL_TRUE, 0, sizeof(vt) * count,
			h_a, 0, NULL, NULL);
	checkError(err, "Copying h_a to device at d_a");

	err = clEnqueueWriteBuffer(commands, d_b, CL_TRUE, 0, sizeof(vt) * count,
			h_b, 0, NULL, NULL);
	checkError(err, "Copying h_b to device at d_b");

	err = clEnqueueWriteBuffer(commands, d_res, CL_TRUE, 0, sizeof(vt),
				&h_res, 0, NULL, NULL);
		checkError(err, "Copying h_b to device at d_b");

	// Set the arguments to our compute kernel
	err = clSetKernelArg(ko_vecdot[0], 0, sizeof(cl_mem), &d_a);
	err |= clSetKernelArg(ko_vecdot[0], 1, sizeof(cl_mem), &d_b);
	err |= clSetKernelArg(ko_vecdot[0], 2, sizeof(cl_mem), &d_c);
	//err |= clSetKernelArg(ko_vadd, 3, sizeof(uint), &npg);
	//err |= clSetKernelArg(ko_vadd, 4, sizeof(uint), &npi);
	err |= clSetKernelArg(ko_vecdot[0], 3, sizeof(uint), &count);
	checkError(err, "Setting kernel arguments");

	err = clSetKernelArg(ko_vecdot[1], 0, sizeof(cl_mem), &d_res);
	err |= clSetKernelArg(ko_vecdot[1], 1, sizeof(cl_mem), &d_c);
	err |= clSetKernelArg(ko_vecdot[1], 2, sizeof(uint), &count);
	checkError(err, "Setting kernel arguments");

	tick_t time = Clock::Tick();
	// Execute the kernel over the entire range of our 1d input data set
	// letting the OpenCL runtime choose the work-group size
	//size_t workgroupsize = 256;
	global = count;
	err = clEnqueueNDRangeKernel(commands, ko_vecdot[0], 1, NULL, &global,
			NULL, 0, NULL, NULL);
	err = clEnqueueNDRangeKernel(commands, ko_vecdot[1], 1, NULL, &global,
			NULL, 0, NULL, NULL);
	//err = clEnqueueNDRangeKernel(commands, ko_vadd, 1, NULL, &global, NULL, 0, NULL, NULL);
	checkError(err, "Enqueueing kernel");

	// Wait for the commands to complete before stopping the timer
	err = clFinish(commands);
	checkError(err, "Waiting for kernel to finish");

	time = Clock::Tick() - time;
	printf("\nThe kernel ran in %lf ms\n", Clock::TicksToMillisecondsF(time));

	// Read back the results from the compute device

	err = clEnqueueReadBuffer(commands, d_res, CL_TRUE, 0, sizeof(vt),
			&h_res, 0, NULL, NULL);
	if (err != CL_SUCCESS) {
		printf("Error: Failed to read output array!\n%s\n", err_code(err));
		exit(1);
	}

	// Test the results
	correct = 0;

	vt sum = 0;
	vt* h_resc = (vt*) calloc(count, sizeof(vt));
	time = Clock::Tick();
	for (i = 0; i < count; i++) {
		h_resc[i] = h_a[i] * h_b[i];     // assign element i of a+b to tmp
		sum += h_resc[i];                // compute deviation of expected and output result
	}
	time = Clock::Tick() - time;
	printf("\nThe CPU ran in %lf ms\n", Clock::TicksToMillisecondsF(time));


	cout << "dot  " << h_res << "  " << sum << endl;

	// summarise results

	// cleanup then shutdown
	clReleaseMemObject(d_a);
	clReleaseMemObject(d_b);
	clReleaseMemObject(d_c);
	clReleaseProgram(program);
	clReleaseKernel(ko_vecdot[0]);
	clReleaseKernel(ko_vecdot[1]);
	clReleaseCommandQueue(commands);
	clReleaseContext(context);

	free(h_a);
	free(h_b);
	free(h_c);

	return 0;
}

}

#endif
