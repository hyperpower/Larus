/************************
 //  \file   CalDef.h
 //  \brief
 // 
 //  \author czhou
 //  \date   14 janv. 2015 
 ***********************/
#ifndef CALDEF_H_
#define CALDEF_H_

#include "../TypeDef.h"
#include "../Grid/SPTree.h"
#include "../Grid/SPTreeNode.h"
#include "../Grid/Forest.h"

namespace Larus {

const LarusDef::size_type Idx_IDX = 0;
const LarusDef::size_type LS_IDX = 0;
const LarusDef::size_type VOF_IDX = 0;
const LarusDef::size_type VOFt_IDX= 0;
const LarusDef::size_type U_IDX   = 2;
const LarusDef::size_type V_IDX   = 3;

template<class NODE>
inline int getIDX(NODE* pn){
	ASSERT(pn!=NULL_PTR);
	return int(pn->data->aCenterData[Idx_IDX]);
}




//work with Forest===============================
//===============================================



//work with tree=================================
//===============================================



}


#endif /* CALDEF_H_ */
