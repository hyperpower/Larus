/************************
 //  \file   Forest.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   6 févr. 2015 
 ***********************/



#include "Forest.h"
#include <iostream>
#include <stdio.h>

namespace Larus {
void is_leaf_has_data(pQTNode pn, utPointer p) {
	if (!pn->hasChild()) {
		if (pn->data != NULL_PTR) {
			(*(CAST(int*, p)))++;
		}
	}
}
void is_leaf_has_data(pOCNode pn, utPointer p) {
	if (!pn->hasChild()) {
		if (pn->data != NULL_PTR) {
			(*(CAST(int*, p)))++;
		}
	}
}

}
