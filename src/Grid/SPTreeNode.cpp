/************************
 //  \file   QuadTreeNode.cpp
 //  \brief
 // 
 //  \author czhou
 //  \date   4 juin 2014 
 ***********************/

#include "SPTreeNode.h"
#include <stdio.h>
#include <iostream>

namespace Larus
{
//====================================================
//SPNode =============================================
SPNodeIdx toSPNodeIdx(int i)
{
	return (i >= 0 && i <= 7) ? SP_NODEIDX[i] : ErrSPIdx;
}

SPNodeIdx toSPNodeIdx(SPDirection i)
{
	return (i >= 0 && i <= 7) ? SP_NODEIDX[i] : ErrSPIdx;
}

const SPNodeBoundary SP_NODEBOU[6] =
{ SPB_W, SPB_N, SPB_E, SPB_S, SPB_F, SPB_B };

SPNodeBoundary toSPNodeBoundary(int i)
{
	return (i >= 4 && i <= 9) ? SP_NODEBOU[int(i) - 4] : ErrSPBoundary;
}

SPNodeBoundary toSPNodeBoundary(SPDirection i)
{
	return (i != ErrSPDirection && i >= 4 && i <= 9) ?
			SP_NODEBOU[int(i) - 4] : ErrSPBoundary;
}
SPDirection toDirection(int i)
{
	return (i >= 0 && i <= 26) ? SP_DIR[i] : ErrSPDirection;
}
SPDirection toDirection(SPNodeIdx i)
{
	return (i != ErrSPIdx) ? SP_DIR[int(i)] : ErrSPDirection;
}
SPDirection toDirection(SPNodeBoundary i)
{
	return (i != ErrSPBoundary) ? SP_DIR[int(i)] : ErrSPDirection;
}

bool isAdjacent(const int &idx_node, const int& idx_boundary)
{
	return SP_ADJ[idx_node][idx_boundary - 4];
}

bool isAdjacent(SPNodeIdx idx, SPNodeBoundary b)
{
	return (idx != ErrSPIdx && b != ErrSPBoundary) ?
			isAdjacent(int(idx), int(b)) : false;
}

int reflecttIdx(int idx_node, int direction)
{
	return SP_REFLECT[idx_node][direction - 4];
}

SPNodeIdx reflecttIdx(SPNodeIdx idx, SPDirection d)
{
	if (d < 4 || d > 9 || idx == ErrSPIdx) {
		std::cerr << "!> SPDirection Index=" << d
				<< " Input error. return false \n";
		std::cerr
				<< "!> Function int reflecttIdx(int idx_node, int direction)\n";
		return ErrSPIdx;
	} else {
		return toSPNodeIdx(reflecttIdx(int(idx), int(d)));
	}
}
SPNodeIdx toQTNodeIdx(int i)
{
	return (i >= 0 && i <= 7) ? SP_NODEIDX[i] : ErrSPIdx;
}
SPDirection reflectDirection(SPDirection d)
{
	return (d != ErrSPDirection) ?
			toDirection(SP_reflectDIR[int(d)]) : ErrSPDirection;
}
SPNodeIdx diagonalIdx(SPNodeIdx idx, SPDirection d)
{
	if (d == ErrSPDirection)
		return ErrSPIdx;
	if (d >= 0 && d <= 3)
		return toSPNodeIdx(SP_DIANODE_XY[int(idx)]);
	if (d >= 10 && d <= 13)
		return toSPNodeIdx(SP_DIANODE_YZ[int(idx)]);
	if (d >= 14 && d <= 17)
		return toSPNodeIdx(SP_DIANODE_ZX[int(idx)]);
	return ErrSPIdx;
}
int diagonalIdx_XY(int idx)
{
	return SP_DIANODE_XY[idx];
}
SPNodeIdx diagonalIdx_XY(SPNodeIdx idx)
{
	return (idx != ErrSPIdx) ? toSPNodeIdx(SP_DIANODE_XY[int(idx)]) : ErrSPIdx;
}

bool dia_Sibling(SPNodeIdx idx, SPDirection d)
{
	if (d == ErrSPDirection)
		return false;
	if (d >= 0 && d <= 3)
		return SP_DIA_SIBLING_XY[int(idx)][int(d)];
	if (d >= 10 && d <= 13)
		return SP_DIA_SIBLING_YZ[int(idx)][int(d) - 10];
	if (d >= 14 && d <= 17)
		return SP_DIA_SIBLING_ZX[int(idx)][int(d) - 14];
	return false;
}
bool dia_Sibling_XY(SPNodeIdx idx, SPDirection d)
{
	//d 0 to 3
	return (d != ErrSPDirection) ? SP_DIA_SIBLING_XY[int(idx)][int(d)] : false;
}
bool dia_Sibling_YZ(SPNodeIdx idx, SPDirection d)
{
	//d 10 to 13
	return (d != ErrSPDirection) ?
			SP_DIA_SIBLING_YZ[int(idx)][int(d) - 10] : false;
}
bool dia_Sibling_ZX(SPNodeIdx idx, SPDirection d)
{
	//d 14 to 17
	return (d != ErrSPDirection) ?
			SP_DIA_SIBLING_ZX[int(idx)][int(d) - 14] : false;
}
bool out_cor(SPNodeIdx idx, SPDirection d)
{
	if (d == ErrSPDirection)
		return false;
	if (d >= 0 && d <= 3)
		return SP_OUT_COR_XY[int(idx)][int(d)];
	if (d >= 10 && d <= 13)
		return SP_OUT_COR_YZ[int(idx)][int(d) - 10];
	if (d >= 14 && d <= 17)
		return SP_OUT_COR_ZX[int(idx)][int(d) - 14];
	return false;
}
bool out_cor_XY(SPNodeIdx idx, SPDirection d)
{
	//d 0 to 3
	return (d != ErrSPDirection) ? SP_OUT_COR_XY[int(idx)][int(d)] : false;
}

bool out_cor_YZ(SPNodeIdx idx, SPDirection d)
{
	//d 10 to 13
	return (d != ErrSPDirection) ? SP_OUT_COR_YZ[int(idx)][int(d) - 10] : false;
}
bool out_cor_ZX(SPNodeIdx idx, SPDirection d)
{
	//d 14 to 17
	return (d != ErrSPDirection) ? SP_OUT_COR_ZX[int(idx)][int(d) - 14] : false;
}
int diagonalIdx_YZ(int idx)
{
	return SP_DIANODE_YZ[idx];
}
SPNodeIdx diagonalIdx_YZ(SPNodeIdx idx)
{
	return (idx != ErrSPIdx) ? toSPNodeIdx(SP_DIANODE_YZ[int(idx)]) : ErrSPIdx;
}
int diagonalIdx_ZX(int idx)
{
	return SP_DIANODE_ZX[idx];
}
SPNodeIdx diagonalIdx_ZX(SPNodeIdx idx)
{
	return (idx != ErrSPIdx) ? toSPNodeIdx(SP_DIANODE_ZX[int(idx)]) : ErrSPIdx;
}
int commonBoundary_XY(int idx1, int idx2)
{
	return SP_COMBOU_XY[idx1][idx2];
}
SPDirection commonDirection(SPNodeIdx idx, SPDirection d)
{
	if (d == ErrSPDirection)
		return ErrSPDirection;
	if (d >= 0 && d <= 3)
		return toDirection(SP_COMDIR_XY[int(idx)][int(d)]);
	if (d >= 10 && d <= 13)
		return toDirection(SP_COMDIR_YZ[int(idx)][int(d) - 10]);
	if (d >= 14 && d <= 17)
		return toDirection(SP_COMDIR_ZX[int(idx)][int(d) - 14]);
	return ErrSPDirection;
}

SPNodeBoundary commonBoundary_XY(SPNodeIdx idx1, SPNodeIdx idx2)
{
	return (idx1 != ErrSPIdx && idx2 != ErrSPIdx) ?
			toSPNodeBoundary(commonBoundary_XY(int(idx1), int(idx2))) :
			ErrSPBoundary;
}
SPDirection commonDirection_XY(SPNodeIdx idx1, SPNodeIdx idx2)
{
	return (idx1 != ErrSPIdx && idx2 != ErrSPIdx) ?
			toDirection(commonBoundary_XY(int(idx1), int(idx2))) :
			ErrSPDirection;
}
SPDirection commonDirection_XY(SPNodeIdx idx1, SPDirection d)
{
	//d 1 to 3
	return (idx1 != ErrSPIdx && d != ErrSPDirection) ?
			toDirection(SP_COMDIR_XY[int(idx1)][int(d)]) : ErrSPDirection;
}
int commonBoundary_YZ(int idx1, int idx2)
{
	return SP_COMBOU_YZ[idx1][idx2];
}
SPNodeBoundary commonBoundary_YZ(SPNodeIdx idx1, SPNodeIdx idx2)
{
	return (idx1 != ErrSPIdx && idx2 != ErrSPIdx) ?
			toSPNodeBoundary(commonBoundary_YZ(int(idx1), int(idx2))) :
			ErrSPBoundary;
}
SPDirection commonDirection_YZ(SPNodeIdx idx1, SPNodeIdx idx2)
{
	return (idx1 != ErrSPIdx && idx2 != ErrSPIdx) ?
			toDirection(commonBoundary_YZ(int(idx1), int(idx2))) :
			ErrSPDirection;
}
SPDirection commonDirection_YZ(SPNodeIdx idx1, SPDirection d)
{
	//d 10 to 13
	return (idx1 != ErrSPIdx && d != ErrSPDirection) ?
			toDirection(SP_COMDIR_YZ[int(idx1)][int(d) - 10]) : ErrSPDirection;
}
int commonBoundary_ZX(int idx1, int idx2)
{
	return SP_COMBOU_ZX[idx1][idx2];
}
SPNodeBoundary commonBoundary_ZX(SPNodeIdx idx1, SPNodeIdx idx2)
{
	return (idx1 != ErrSPIdx && idx2 != ErrSPIdx) ?
			toSPNodeBoundary(commonBoundary_ZX(int(idx1), int(idx2))) :
			ErrSPBoundary;
}
SPDirection commonDirection_ZX(SPNodeIdx idx1, SPNodeIdx idx2)
{
	return (idx1 != ErrSPIdx && idx2 != ErrSPIdx) ?
			toDirection(commonBoundary_ZX(int(idx1), int(idx2))) :
			ErrSPDirection;
}
SPDirection commonDirection_ZX(SPNodeIdx idx1, SPDirection d)
{
	//d 14 to 17
	return (idx1 != ErrSPIdx && d != ErrSPDirection) ?
			toDirection(SP_COMDIR_ZX[int(idx1)][int(d) - 14]) : ErrSPDirection;
}

SPDirection transfer_3Ddir_to_2Ddir(SPNodeIdx idx, SPDirection d)
{
	//d 18 to 25
	return (idx != ErrSPIdx && d != ErrSPDirection) ?
			toDirection(SP_DIR_3Dto2D[int(idx)][int(d) - 18]) : ErrSPDirection;
}
int Direction_Decompose(const SPDirection& ind, SPDirection& d1,
		SPDirection& d2, SPDirection& d3)
{
	d1 = toDirection(SP_DIR_DECOMPOSE[ind][1]);
	d2 = toDirection(SP_DIR_DECOMPOSE[ind][2]);
	d3 = toDirection(SP_DIR_DECOMPOSE[ind][3]);
	return SP_DIR_DECOMPOSE[ind][0];
}

SPDirection Direction_Compose(const SPDirection& dx, const SPDirection& dy)
{
	//dx 4 or 6
	//dy 5 or 7
	ASSERT(dx != ErrSPDirection && dy != ErrSPDirection);
	if (dx == SPD_IM) {
		return (dy == SPD_JM) ? SPD_MM : SPD_MP;
	} else { //dx==SPD_IP
		return (dy == SPD_JM) ? SPD_PM : SPD_PP;
	}
}
//SPNode ----------------------------------------

//OCNode=========================================

//QTNode ==============================================

void draw_boundary(pQTNode qnode, FILE* file)
{
	fprintf(file, "%f %f\n", qnode->cell->getMM().x, qnode->cell->getMM().y);
	fprintf(file, "%f %f\n\n", qnode->cell->getMP().x, qnode->cell->getMP().y);
	fprintf(file, "%f %f\n", qnode->cell->getMP().x, qnode->cell->getMP().y);
	fprintf(file, "%f %f\n\n", qnode->cell->getPP().x, qnode->cell->getPP().y);
	fprintf(file, "%f %f\n", qnode->cell->getPP().x, qnode->cell->getPP().y);
	fprintf(file, "%f %f\n\n", qnode->cell->getPM().x, qnode->cell->getPM().y);
	fprintf(file, "%f %f\n", qnode->cell->getPM().x, qnode->cell->getPM().y);
	fprintf(file, "%f %f\n\n", qnode->cell->getMM().x, qnode->cell->getMM().y);
}

//QTNode ----------------------------------------------
pQTNode getChild_vertex(pQTNode pn, SPDirection dir)
{
	ASSERT(0 <= int(dir) && int(dir) <= 3);
	if (pn->child[int(dir)] != NULL_PTR) {
		pn = pn->child[int(dir)];
		return getChild_vertex(pn, dir);
	} else if (pn->hasChild()) {
		return NULL_PTR;
	} else {
		return pn;
	}
}
//=====================================================
bool hasType(int type, SPNodeType t)
{
	return (type | t) == type ? true : false;
}

std::string parseQTNodeType(int type)
{
	int flag = 0;
	if (hasType(type, SPT_ghost)) {
		return " ghost ";
	}
	if (hasType(type, SPT_normal)) {
		return " normal ";
	}
	if (hasType(type, SPT_cut)) {
		return " cut ";
	}
	if (flag == 0) {
		return " Error input node type";
	}
	return " Error input node type";
}

std::string parseSPNodeFaceType(const SPNodeFaceType& t)
{
	switch (t) {
	case SPFT_Boundary:
		return "SPFT_Boundary";
		break;
	case SPFT_Equal:
		return "SPFT_Equal";
		break;
	case SPFT_FineCoarse:
		return "SPFT_FineCorase";
		break;
	case SPFT_CoarseFine:
		return "SPFT_CoraseFine";
		break;
	default:
		break;
	}
	return "SPFT_Error";
}

int point_at_which_child(pQTNode pnode, Point2D& p)
{
	if (pnode != NULL) {
		return pnode->whichChild(p);
	}
	return -1;
}
int point_at_which_child(pOCNode pnode, Point3D& p)
{
	if (pnode != NULL) {
		return pnode->whichChild(p);
	}
	return -1;
}
template<typename NODE>
int T_relation_node_level(NODE* pnode1, NODE* pnode2)
{
	if (pnode1 == NULL_PTR && pnode1 == NULL_PTR) {
		return -3;
	}
	if (pnode1 == NULL_PTR) {
		return -1;
	}
	if (pnode1 == NULL_PTR) {
		return 1;
	}
	if (pnode1->getLevel() == pnode2->getLevel()) {
		return 0;
	}
	if (pnode1->getLevel() > pnode2->getLevel()) {
		return -2;
	} else {
		return 2;
	}
}

int relation_node_level(pQTNode pnode1, pQTNode pnode2)
{
	return T_relation_node_level(pnode1, pnode2);
}
int relation_node_level(pOCNode pnode1, pOCNode pnode2)
{
	return T_relation_node_level(pnode1, pnode2);
}

} //this is the end of namespace
