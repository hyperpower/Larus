/************************
 //  \file   Levelset.h
 //  \brief
 // 
 //  \author czhou
 //  \date   10 févr. 2015 
 ***********************/
#ifndef LEVELSET_H_
#define LEVELSET_H_

#include "../TypeDef.h"
#include "CalDef.h"
#include "Levelset.h"
#include "../Grid/SPTree.h"
#include "../Grid/SPTreeNode.h"

namespace Larus {
//work with tree=================================
//===============================================
//set value======================================
void set_sphere_ls(pOCTree oct, Point3D p, Float r, LarusDef::size_type idx=LS_IDX);

void set_sphere_ls_initial_refine(pOCTree oct,
		Point3D p,
		Float r,
		LarusDef::size_type idx=LS_IDX);
}



#endif /* LEVELSET_H_ */
