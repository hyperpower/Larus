#pragma OPENCL EXTENSION cl_khr_fp64 : enable
__kernel void vec_add(__global double *A, 
  __global double *B,
  __global double *C,
   unsigned int count){                                                                                      
   // Get the work-item’s unique ID             
   uint idx = get_global_id(0);                  
                                                
   // Add the corresponding locations of        
   // 'A' and 'B', and store the result in 'C'. 
   if(idx < count){
        C[idx] = A[idx] + B[idx];
   }                    
}          

__kernel void vec_times_num(__global const double *num,
                            __global const double *vec,
                            __global double *vec_res,
                            unsigned int count){                                                                                      
   // Get the work-item’s unique ID             
   uint idx = get_global_id(0);                  
                                                
   // Add the corresponding locations of        
   // 'A' and 'B', and store the result in 'C'. 
   if(idx < count){
        vec_res[idx] = num[0] * vec[idx];
   }                    
}

__kernel void vec_sum(__global double *sum,
                      __global const double *vec,
                      uint count){                                                                                      
   // Get the work-item’s unique ID             
   int idx = get_global_id(0);                  
                                                
   // Add the corresponding locations of        
   // 'A' and 'B', and store the result in 'C'. 
   if(idx < count){
        sum[0] += vec[idx];
   }                    
}

__kernel void vec_dot(
      __global const double* x, // input vector
      __global const double* y, // input vector
      __global double *r, // result vector
      int n // input vector size
){
    int id = get_global_id(0);
    if ( id < n ){
        r[id] = x[id] * y[id]; // multiply elements, store product
    }
}

#define LOCAL_GROUP_XDIM 256
// Kernel for part 1 of dot product, version 3.
__kernel __attribute__((reqd_work_group_size(LOCAL_GROUP_XDIM, 1, 1)))
void dot_persist_kernel(
    __global const double * x, // input vector
    __global const double * y, // input vector
    __global double * r,       // result vector
    uint n_per_group,          // elements processed per group
    uint n_per_work_item,      // elements processed per work item
    uint n                     // input vector size
    ){
    // uint id = get_global_id(0); // unused here
    uint lcl_id = get_local_id(0);
    uint grp_id = get_group_id(0);
    double priv_acc = 0; // accumulator in private memory
    __local double lcl_acc[LOCAL_GROUP_XDIM]; // accumulators in local memory
    uint grp_off = mul24(n_per_group, grp_id); // group offset
    uint lcl_off = grp_off + lcl_id; // local offset

    // Accumulate products over n_per_work_item elements.
    double priv_val = 0;
    for ( uint i = 0; i < n_per_work_item; i++, lcl_off += LOCAL_GROUP_XDIM)
    {
        // Be wary of out of range offsets, just add 0 if out of range.
        // This code uses conditional expressions rather than ifs for efficiency.
        bool in_range = ( lcl_off < n );
        uint lcl_off2 = ( in_range ) ? lcl_off : 0;
        priv_val = x[lcl_off2] * y[lcl_off2]; // multiply elements
        priv_acc += ( in_range ) ? priv_val : 0; // accumulate result
    }

    // Store result accumulated so far to local accumulator.
    lcl_acc[lcl_id] = priv_acc;
    barrier(CLK_LOCAL_MEM_FENCE);
    // Find the sum of the accumulated products.
    uint dist = LOCAL_GROUP_XDIM; // i.e., get_local_size(0);
    while ( dist > 1 )
    {
        dist >>= 1;
        if ( lcl_id < dist ){
            // Private memory accumulator avoids extra local memory read.
            priv_acc += lcl_acc[lcl_id + dist];
            lcl_acc[lcl_id] = priv_acc;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // Store the result.
    if ( lcl_id == 0 ){
        r[grp_id] = priv_acc;
    }

}


// Kernel for part 1 of dot product, version 2, using local reduction.
__kernel __attribute__((reqd_work_group_size(LOCAL_GROUP_XDIM, 1, 1)))
void dot_local_reduce_kernel(
      __global const double * x, // input vector
      __global const double * y, // input vector
      __global double * r, // result vector
      uint n // input vector size
){
  uint id     = get_global_id(0);
  uint lcl_id = get_local_id(0);
  uint grp_id = get_group_id(0);
  double priv_acc = 0;                      // accumulator in private memory
  __local double lcl_acc[LOCAL_GROUP_XDIM]; // accumulators in local memory
  if ( id < n ){
      priv_acc = lcl_acc[lcl_id] = x[id] * y[id]; // multiply elements, store product
  }
  barrier(CLK_LOCAL_MEM_FENCE); // Find the sum of the accumulators.
  uint dist = LOCAL_GROUP_XDIM; // i.e., get_local_size(0);
  while ( dist >= 1 ){
      dist >>= 1;
      if ( lcl_id < dist ){
          // Private memory accumulator avoids extra local memory read.
          priv_acc += lcl_acc[lcl_id + dist];
          lcl_acc[lcl_id] = priv_acc;
      }
      barrier(CLK_LOCAL_MEM_FENCE);
  }
  // Store the result (the sum for the local work group).
  if ( lcl_id == 0 ){
      r[grp_id] = priv_acc;
  }
}



