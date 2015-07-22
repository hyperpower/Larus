/************************
 //  \file   TypeDef.cpp
 //  \brief
 // 
 //  \author zhou
 //  \date   23 mai 2014 
 ***********************/


#include "TypeDef.h"
#include <sys/time.h>
#include <time.h>

namespace Larus{
bool isEqual(Float a,Float b){
	if(ABS(a-b)<=SMALL){
		return true;
	}else{
		return false;
	}
}

bool isZero(Float a){
	if(ABS(a)<SMALL){
		return true;
	}else{
		return false;
	}
}

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

}
