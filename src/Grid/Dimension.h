/************************
 //  \file   Dimension.h
 //  \brief
 // 
 //  \author czhou
 //  \date   5 févr. 2015 
 ***********************/
#ifndef DIMENSION_H_
#define DIMENSION_H_

#include "../TypeDef.h"

#include "SPTree.h"
#include "SPTreeNode.h"
#include "CellData.h"
#include "../Algebra/Matrix.h"
#include "../Algebra/Space.h"

namespace Larus {

class Dimension_3D {
public:
	typedef LarusDef::size_type size_type;
	static const size_type DIM = 3;
	static const size_type NUM_CELLS = 8;
	static const size_type NUM_NEBOR = 6;
	typedef SPTree<OCNode, 3> Tree;
	typedef SPTree<OCNode, 3>* pTree;
	typedef OCNode Node;
	typedef OCNode* pNode;
	typedef Point3D Point;
	typedef Point3D* pPoint;
	typedef Cell3D Cell;
	typedef Cell3D* pCell;
	typedef CellData3D CellData;
	typedef CellData3D* pCellData;
	typedef CellData::value_type value_type;
};

class Dimension_2D {
public:
	typedef LarusDef::size_type size_type;
	static const size_type DIM = 2;
	static const size_type NUM_CELLS = 4;
	static const size_type NUM_NEBOR = 4;
	typedef SPTree<QTNode, 2> Tree;
	typedef SPTree<QTNode, 2>* pTree;
	typedef QTNode  Node;
	typedef QTNode* pNode;
	typedef Point2D Point;
	typedef Point2D* pPoint;
	typedef Cell2D Cell;
	typedef Cell2D* pCell;
	typedef CellData2D CellData;
	typedef CellData2D* pCellData;
	typedef CellData::value_type value_type;
};

}

#endif /* DIMENSION_H_ */
