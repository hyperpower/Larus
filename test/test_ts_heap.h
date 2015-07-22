/*
 * test_ts_heap.h
 *
 *  Created on: May 30, 2015
 *      Author: zhou
 */

#ifndef TEST_TS_HEAP_H_
#define TEST_TS_HEAP_H_

#include "../src/TS/ts_vertex.h"
#include "../src/TS/ts_point.h"
namespace LarusTS{

bool compare(int a , int b){
	return a>b;
}

void test_test_ts(){
    Point<double, 3> a(0,0,0);
    Point<double, 3> b(1,1,1);
    Point<double, 3> c(3,3,2);
    Point<double, 3> d(0,0,1);
    Vertex<double, 3> ver(0,1,2);
    ver.show();

    cout<< point_orientation(a,b,c)<<endl;
    cout<< point_orientation_sos(a,b,c)<<endl;
}
}



#endif /* TEST_TS_HEAP_H_ */
