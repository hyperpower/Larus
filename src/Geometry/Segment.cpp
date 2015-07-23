/************************
 //  \file   Segment.cpp
 //  \brief
 // 
 //  \author zhou
 //  \date   23 mai 2014 
 ***********************/
#include "Segment.h"
#include "Point.h"
#include <assert.h>
#include <math.h>
#include <list>
#include <iostream>

namespace Larus {

void parseIntersectType(int inter) {
	switch (inter) {
	case NO_INTERSECT:
		std::cout << "NO_INTERSERSECT \n";
		break;
	case INTERSECT:
		std::cout << "INTERSERSECT normal\n";
		break;
	case INTERSECT | START_1:
		std::cout << "INTERSERSECT 1 with 1_START\n";
		break;
	case INTERSECT | START_2:
		std::cout << "INTERSERSECT 1 with 2_START\n";
		break;
	case INTERSECT | END_1:
		std::cout << "INTERSERSECT 1 with 1_END\n";
		break;
	case INTERSECT | END_2:
		std::cout << "INTERSERSECT 1 with 2_END\n";
		break;
	case INTERSECT | START_1 | START_2:
		std::cout << "INTERSERSECT 2 : 1_2_S\n";
		break;
	case INTERSECT | START_1 | END_2:
		std::cout << "INTERSERSECT 2 : 1_S 2_E\n";
		break;
	case INTERSECT | START_2 | END_1:
		std::cout << "INTERSERSECT 2 : 1_E 2_S\n";
		break;
	case INTERSECT | END_2 | END_1:
		std::cout << "INTERSERSECT 2 : 1_2_E \n";
		break;
	default:
		std::cout << "ERR parse\n";
		break;
	}
}

//this is the end of the class Segment2D
//----------------------------
    

}
//end namespace

