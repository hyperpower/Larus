/************************
 //  \file   QuadTreeNode.h
 //  \brief
 // 
 //  \author czhou
 //  \date   4 juin 2014 
 ***********************/
#ifndef QUADTREENODE_H_
#define QUADTREENODE_H_

#include "../TypeDef.h"
#include "../Utility/ArrayList.h"
#include "CellData.h"
#include "Cell.h"
#include "../IO/IO.h"
#include "../IO/IO_vtk.h"

namespace Larus {

//===============================================
//Start of SPNode================================
//===============================================

enum SPNodeType {
	SPT_normal = 1 << 0, SPT_ghost = 1 << 1, SPT_cut = 1 << 2,
};
enum SPNodeIdx {
	//=========================
	//
	//   ---------------
	//   |  MP  |  PP  |
	//   |  0   |  2   |
	//   |  WN  |  NE  |
	//   ---------------
	//   |  SW  |  SE  |
	//   |  1   |  3   |
	//   |  MM  |  PM  |
	//   ---------------
	//
	//   ---------------  ---------------
	//   |  MPM |  PPM |  |  MPP |  PPP |
	//   |  0   |  2   |  |  4   |  6   |
	//   |  WNB |  NEB |  |  WNF |  NEF |
	//   ---------------  ---------------
	//   |  SWB |  SEB |  |  SWP |  SEP |
	//   |  1   |  3   |  |  5   |  7   |
	//   |  MMM |  PMM |  |  MMM |  PMM |
	//   ---------------  ---------------
	//=========================
	ErrSPIdx = -1,
	//2D=======================
	PP2 = 2,
	MM2 = 1,
	PM2 = 3,
	MP2 = 0,
	NW2 = 0,
	SW2 = 1,
	NE2 = 2,
	SE2 = 3,
	//3D=======================
	PPM = 2,
	MMM = 1,
	PMM = 3,
	MPM = 0,
	NWB = 0,
	SWB = 1,
	NEB = 2,
	SEB = 3,
	//3D=====Front=============
	PPP = 6,
	MMP = 5,
	PMP = 7,
	MPP = 4,
	NWF = 4,
	SWF = 5,
	NEF = 6,
	SEF = 7,
};

enum SPDirection {
	ErrSPDirection = -3,
	//2D=======================
	SPD_IP = 6,
	SPD_IM = 4,
	SPD_JP = 5,
	SPD_JM = 7,

	SPD_MP = 0,
	SPD_MM = 1,
	SPD_PP = 2,
	SPD_PM = 3,

	SPD_E = 6,
	SPD_W = 4,
	SPD_N = 5,
	SPD_S = 7,

	SPD_NE = 2,
	SPD_SW = 1,
	SPD_SE = 3,
	SPD_NW = 0,
	//W -> N -> E -> S
	//4    5    6    7

	//3D=======================
	SPD_KP = 8,
	SPD_KM = 9,

	SPD_F = 8,
	SPD_B = 9,

	SPD_PP_XY = 2,
	SPD_MM_XY = 1,
	SPD_PM_XY = 3,
	SPD_MP_XY = 0,

	SPD_PP_YZ = 2 + 10,
	SPD_MM_YZ = 1 + 10,
	SPD_PM_YZ = 3 + 10,
	SPD_MP_YZ = 0 + 10,

	SPD_PP_ZX = 2 + 14,
	SPD_MM_ZX = 1 + 14,
	SPD_PM_ZX = 3 + 14,
	SPD_MP_ZX = 0 + 14,
	//SPD_XYZ==================
	SPD_MMM = 17 + 1,
	SPD_PMM = 17 + 2,
	SPD_MPM = 17 + 3,
	SPD_PPM = 17 + 4,
	SPD_MMP = 17 + 5,
	SPD_PMP = 17 + 6,
	SPD_MPP = 17 + 7,
	SPD_PPP = 17 + 8,
};

enum SPPlane {
	ErrSPPlane = -4, SP_XY = 0, SP_YZ = 10, SP_ZX = 14,
};

const SPDirection SP_DIR[26] = { SPD_MP_XY, SPD_MM_XY, SPD_PP_XY, SPD_PM_XY,
		SPD_W, SPD_N, SPD_E, SPD_S, SPD_KP, SPD_KM, SPD_MP_YZ, SPD_MM_YZ,
		SPD_PP_YZ, SPD_PM_YZ, SPD_MP_ZX, SPD_MM_ZX, SPD_PP_ZX, SPD_PM_ZX,
		SPD_MMM, SPD_PMM, SPD_MPM, SPD_PPM, SPD_MMP, SPD_PMP, SPD_MPP, SPD_PPP };

enum SPNodeBoundary {
	ErrSPBoundary = -2,
	SPB_IP = 6,
	SPB_IM = 4,
	SPB_JP = 5,
	SPB_JM = 7,
	SPB_KP = 8,
	SPB_KM = 9,
	SPB_E = 6,
	SPB_W = 4,
	SPB_N = 5,
	SPB_S = 7,
	SPB_F = 8,
	SPB_B = 9,

};
const SPDirection SP_DIR_OP[26] = //
		{ //
		SPD_PM_XY,  //SPD_MP = 0,
				SPD_PP_XY,  //SPD_MM = 1,
				SPD_MM_XY,  //SPD_PP = 2,
				SPD_MP_XY,  //SPD_PM = 3,
				SPD_IP,     //SPD_IM = 4,
				SPD_JM,     //SPD_JP = 5,
				SPD_IM,     //SPD_IP = 6,
				SPD_JP,     //SPD_JM = 7,
				SPD_KM,     //SPD_KP = 8,
				SPD_KP,     //SPD_KM = 9,
				SPD_PM_YZ,  //SPD_MP = 10 + 0,
				SPD_PP_YZ,  //SPD_MM = 10 + 1
				SPD_MM_YZ,  //SPD_PP = 10 + 2,
				SPD_MP_YZ,  //SPD_PM = 10 + 3,
				SPD_PM_ZX,  //SPD_MP = 14 + 0,
				SPD_PP_ZX,  //SPD_MM = 14 + 1,
				SPD_MM_ZX,  //SPD_PP = 14 + 2,
				SPD_MP_ZX,  //SPD_PM = 14 + 3,
				SPD_PPP,    //SPD_MMM = 17 + 1,
				SPD_MPP,    //SPD_PMM = 17 + 2,
				SPD_PMP,    //SPD_MPM = 17 + 3,
				SPD_MMP,    //SPD_PPM = 17 + 4,
				SPD_PPM,    //SPD_MMP = 17 + 5,
				SPD_MPM,    //SPD_PMP = 17 + 6,
				SPD_PMM,    //SPD_MPP = 17 + 7,
				SPD_MMM,    //SPD_PPP = 17 + 8,
		};
const int SP_DIR_MorP[26] =    //
		{    //
		        -1,//SPD_MP = 0,
				-1,//SPD_MM = 1,
				1,//SPD_PP = 2,
				-1,//SPD_PM = 3,
				-1,//SPD_IM = 4,
				1,//SPD_JP = 5,
				1,//SPD_IP = 6,
				-1,//SPD_JM = 7,
				1,//SPD_KP = 8,
				-1,//SPD_KM = 9,
				-1,//SPD_MP = 10 + 0,
				-1,//SPD_MM = 10 + 1
				1,//SPD_PP = 10 + 2,
				-1,//SPD_PM = 10 + 3,
				-1,//SPD_MP = 14 + 0,
				-1,//SPD_MM = 14 + 1,
				1,//SPD_PP = 14 + 2,
				-1,//SPD_PM = 14 + 3,
				-1,//SPD_MMM = 17 + 1,
				-1,//SPD_PMM = 17 + 2,
				-1,//SPD_MPM = 17 + 3,
				-1,//SPD_PPM = 17 + 4,
				-1,//SPD_MMP = 17 + 5,
				-1,//SPD_PMP = 17 + 6,
				-1,//SPD_MPP = 17 + 7,
				1,//SPD_PPP = 17 + 8,
		};

const SPNodeIdx SP_NODEIDX[8] = { MPM, MMM, PPM, PMM, MPP, MMP, PPP, PMP };

const bool SP_ADJ[8][6] = { //
		{ true, true, false, false, false, true }, //
				{ true, false, false, true, false, true }, //
				{ false, true, true, false, false, true }, //
				{ false, false, true, true, false, true }, //
				{ true, true, false, false, true, false }, //
				{ true, false, false, true, true, false }, //
				{ false, true, true, false, true, false }, //
				{ false, false, true, true, true, false } };

const int SP_REFLECT[8][6] = { //
							   //
		{ 2, 1, 2, 1, 4, 4 }, //
				{ 3, 0, 3, 0, 5, 5 }, //
				{ 0, 3, 0, 3, 6, 6 }, //
				{ 1, 2, 1, 2, 7, 7 }, //
				{ 6, 5, 6, 5, 0, 0 }, //
				{ 7, 4, 7, 4, 1, 1 }, //
				{ 4, 7, 4, 7, 2, 2 }, //
				{ 5, 6, 5, 6, 3, 3 } };

const int SP_reflectDIR[26] = { 3, 2, 1, 0, 6, 7, 4, 5, 9, 8, 13, 12, 11, 10,
		17, 16, 15, 14, 25, 24, 23, 22, 21, 20, 19, 18 };

const int SP_DIANODE_XY[8] = { 3, 2, 1, 0, 7, 6, 5, 4 };
const int SP_DIANODE_YZ[8] = { 5, 4, 7, 6, 1, 0, 3, 2 };
const int SP_DIANODE_ZX[8] = { 6, 7, 4, 5, 2, 3, 0, 1 };
const int SP_DIANODE_XYZ[8] = { 7, 6, 5, 4, 3, 2, 1, 0 };

const bool SP_DIA_SIBLING_XY[8][4] = //
		{ { false, false, false, true }, //
				{ false, false, true, false }, //
				{ false, true, false, false }, //
				{ true, false, false, false }, //
				{ false, false, false, true }, //
				{ false, false, true, false }, //
				{ false, true, false, false }, //
				{ true, false, false, false }, //
		};
const bool SP_DIA_SIBLING_YZ[8][4] = //
		{ { true, false, false, false }, //
				{ false, false, true, false }, //
				{ true, false, false, false }, //
				{ false, false, true, false }, //
				{ false, true, false, false }, //
				{ false, false, false, true }, //
				{ false, true, false, false }, //
				{ false, false, false, true }  //
		};
const bool SP_DIA_SIBLING_ZX[8][4] = //
		{ { false, false, true, false }, //
				{ false, false, true, false }, //
				{ false, false, false, true }, //
				{ false, false, false, true }, //
				{ true, false, false, false }, //
				{ true, false, false, false }, //
				{ false, true, false, false }, //
				{ false, true, false, false }  //
		};
const bool SP_OUT_COR_XY[8][4] = {
//
		{ true, false, false, false }, //
		{ false, true, false, false }, //
		{ false, false, true, false }, //
		{ false, false, false, true }, //
		{ true, false, false, false }, //
		{ false, true, false, false }, //
		{ false, false, true, false }, //
		{ false, false, false, true } //
};

const bool SP_OUT_COR_YZ[8][4] = //
		{ { false, false, false, true }, //
				{ false, true, false, false }, //
				{ false, false, false, true }, //
				{ false, true, false, false }, //
				{ false, false, true, false }, //
				{ true, false, false, false }, //
				{ false, false, true, false }, //
				{ true, false, false, false } //
		};
const bool SP_OUT_COR_ZX[8][4] = //
		{ { false, true, false, false }, //
				{ false, true, false, false }, //
				{ true, false, false, false }, //
				{ true, false, false, false }, //
				{ false, false, false, true }, //
				{ false, false, false, true }, //
				{ false, false, true, false }, //
				{ false, false, true, false } //
		};
const int SP_COMDIR_XY[8][4] = { //
		{ -2, 4, 5, -2 }, //
				{ 4, -2, -2, 7 }, //
				{ 5, -2, -2, 6 }, //
				{ -2, 7, 6, -2 }, //
				{ -2, 4, 5, -2 }, //
				{ 4, -2, -2, 7 }, //
				{ 5, -2, -2, 6 }, //
				{ -2, 7, 6, -2 } };
const int SP_COMDIR_YZ[8][4] = { //
		{ -2, 9, 5, -2 }, //
				{ 7, -2, -2, 9 }, //
				{ -2, 9, 5, -2 }, //
				{ 7, -2, -2, 9 }, //
				{ 8, -2, -2, 5 }, //
				{ -2, 7, 8, -2 }, //
				{ 8, -2, -2, 5 }, //
				{ -2, 7, 8, -2 } };

const int SP_COMDIR_ZX[8][4] = { //
		{ 9, -2, -2, 4 }, //
				{ 9, -2, -2, 4 }, //
				{ -2, 9, 6, -2 }, //
				{ -2, 9, 6, -2 }, //
				{ -2, 4, 8, -2 }, //
				{ -2, 4, 8, -2 }, //
				{ 6, -2, -2, 8 }, //
				{ 6, -2, -2, 8 } };

const int SP_COMBOU_XY[8][8] = //
		{ { -2, 4, 5, -2, -2, -2, -2, -2 }, //
				{ 4, -2, -2, 7, -2, -2, -2, -2 }, //
				{ 5, -2, -2, 6, -2, -2, -2, -2 }, //
				{ -2, 7, 6, -2, -2, -2, -2, -2 }, //
				{ -2, -2, -2, -2, -2, 4, 5, -2 }, //
				{ -2, -2, -2, -2, 4, -2, -2, 7 }, //
				{ -2, -2, -2, -2, 5, -2, -2, 6 }, //
				{ -2, -2, -2, -2, -2, 7, 6, -2 } };
const int SP_COMBOU_YZ[8][8] = //
		{ { -2, 9, -2, -2, 5, -2, -2, -2 }, //
				{ 9, -2, -2, -2, -2, 7, -2, -2 }, //
				{ -2, -2, -2, 9, -2, -2, 5, -2 }, //
				{ -2, -2, 9, -2, -2, -2, -2, 7 }, //
				{ 5, -2, -2, -2, -2, 8, -2, -2 }, //
				{ -2, 7, -2, -2, 8, -2, -2, -2 }, //
				{ -2, -2, 5, -2, -2, -2, -2, 8 }, //
				{ -2, -2, -2, 7, -2, -2, 8, -2 } };
const int SP_COMBOU_ZX[8][8] = //
		{ { -2, -2, 9, -2, 4, -2, -2, -2 }, //
				{ -2, -2, -2, 9, -2, 4, -2, -2 }, //
				{ 9, -2, -2, -2, -2, -2, 6, -2 }, //
				{ -2, 9, -2, -2, -2, -2, -2, 6 }, //
				{ 4, -2, -2, -2, -2, -2, 8, -2 }, //
				{ -2, 4, -2, -2, -2, -2, -2, 8 }, //
				{ -2, -2, 6, -2, 8, -2, -2, -2 }, //
				{ -2, -2, -2, 6, -2, 8, -2, -2 } };
const bool SP_VERTEX_XYZ[8][8] = {
//          0     1       2        3      4       5      6       7
		{ false, false, true, false, false, false, false, false }, // 0
		{ true, false, false, false, false, false, false, false }, // 1
		{ false, false, false, true, false, false, false, false }, // 2
		{ false, true, false, false, false, false, false, false }, // 3
		{ false, false, false, false, false, false, true, false }, // 4
		{ false, false, false, false, true, false, false, false }, // 5
		{ false, false, false, false, false, false, false, true }, // 6
		{ false, false, false, false, false, true, false, false } // 7
};
const int SP_DIR_3Dto2D[8][8] = {
//         0    1  2  3   4   5   6   7
		{ SP_ZX + 1, 9, 2, SP_YZ + 3, 4, -2, SP_XY + 0, 5 }, //
		{ 0, SP_YZ + 1, SP_ZX + 1, 9, SP_XY + 1, 7, 4, -2 }, //
		{ 9, SP_ZX + 0, SP_YZ + 3, 3, -2, 6, 5, SP_XY + 2 }, //
		{ SP_YZ + 1, 1, 9, SP_ZX + 0, 7, SP_XY + 3, -2, 6 }, //
		{ 4, -2, SP_XY + 0, 5, SP_ZX + 3, 8, 6, SP_YZ + 2 }, //
		{ SP_XY + 1, 7, 4, -2, 4, SP_YZ + 0, SP_ZX + 3, 8 }, //
		{ -2, 6, 5, SP_XY + 2, 8, SP_ZX + 2, SP_YZ + 2, 7 }, //
		{ 7, SP_XY + 3, -2, 6, SP_YZ + 0, 5, 8, SP_ZX + 2 }, //
		};
const int SP_DIR_DECOMPOSE[26][4] = { //
									  //
		{ 2, SPD_IM, SPD_JP, ErrSPDirection },	//SPD_MP = 0,
				{ 2, SPD_IM, SPD_JM, ErrSPDirection },	//SPD_MM = 1,
				{ 2, SPD_IP, SPD_JP, ErrSPDirection },	//SPD_PP = 2,
				{ 2, SPD_IP, SPD_JM, ErrSPDirection },	//SPD_PM = 3,
				{ 1, SPD_IM, ErrSPDirection, ErrSPDirection },	//SPD_IM = 4,
				{ 1, SPD_JP, ErrSPDirection, ErrSPDirection },	//SPD_JP = 5,
				{ 1, SPD_IP, ErrSPDirection, ErrSPDirection },	//SPD_IP = 6,
				{ 1, SPD_JM, ErrSPDirection, ErrSPDirection },	//SPD_JM = 7,
				{ 1, SPD_KP, ErrSPDirection, ErrSPDirection },	//SPD_KP = 8,
				{ 1, SPD_KM, ErrSPDirection, ErrSPDirection },	//SPD_KM = 9,
				{ 2, SPD_JM, SPD_KP, ErrSPDirection },	//SPD_MP = 10 + 0,
				{ 2, SPD_JM, SPD_KM, ErrSPDirection },	//SPD_MM = 10 + 1
				{ 2, SPD_JP, SPD_KP, ErrSPDirection },	//SPD_PP = 10 + 2,
				{ 2, SPD_JP, SPD_KM, ErrSPDirection },	//SPD_PM = 10 + 3,
				{ 2, SPD_KM, SPD_IP, ErrSPDirection },	//SPD_MP = 14 + 0,
				{ 2, SPD_KM, SPD_IM, ErrSPDirection },	//SPD_MM = 14 + 1,
				{ 2, SPD_KP, SPD_IP, ErrSPDirection },	//SPD_PP = 14 + 2,
				{ 2, SPD_KP, SPD_IM, ErrSPDirection },	//SPD_PM = 14 + 3,
				{ 3, SPD_IM, SPD_JM, SPD_KM },		       //SPD_MMM = 17 + 1,
				{ 3, SPD_IP, SPD_JM, SPD_KM },		       //SPD_PMM = 17 + 2,
				{ 3, SPD_IM, SPD_JP, SPD_KM },		       //SPD_MPM = 17 + 3,
				{ 3, SPD_IP, SPD_JP, SPD_KM },		       //SPD_PPM = 17 + 4,
				{ 3, SPD_IM, SPD_JM, SPD_KP },		       //SPD_MMP = 17 + 5,
				{ 3, SPD_IP, SPD_JM, SPD_KP },		       //SPD_PMP = 17 + 6,
				{ 3, SPD_IM, SPD_JP, SPD_KP },		       //SPD_MPP = 17 + 7,
				{ 3, SPD_IP, SPD_JP, SPD_KP }		       //SPD_PPP = 17 + 8,
		};

//Function=======================================
// Bool =========================================
bool isAdjacent(const int &idx_node, const int& idx_boundary);
bool isAdjacent(SPNodeIdx idx, SPNodeBoundary b);
inline bool isFaceDirection(const SPDirection& d, const int& dim) {
	ASSERT(dim == 3 || dim == 2);
	return (dim == 2) ? (d >= 4 && d <= 7) : (d >= 4 && d <= 9);
}

int reflecttIdx(int idx_node, int direction);
SPNodeIdx reflecttIdx(SPNodeIdx idx, SPDirection d);

SPDirection reflectDirection(SPDirection);

SPNodeIdx diagonalIdx(SPNodeIdx idx, SPDirection d);
int diagonalIdx_XY(int idx);
SPNodeIdx diagonalIdx_XY(SPNodeIdx idx);
int diagonalIdx_YZ(int idx);
SPNodeIdx diagonalIdx_YZ(SPNodeIdx idx);
int diagonalIdx_ZX(int idx);
SPNodeIdx diagonalIdx_ZX(SPNodeIdx idx);

bool dia_Sibling(SPNodeIdx idx, SPDirection d);
bool dia_Sibling_XY(SPNodeIdx idx, SPDirection d);
bool dia_Sibling_YZ(SPNodeIdx idx, SPDirection d);
bool dia_Sibling_ZX(SPNodeIdx idx, SPDirection d);

bool out_cor(SPNodeIdx idx, SPDirection d);
bool out_cor_XY(SPNodeIdx idx, SPDirection d);
bool out_cor_YZ(SPNodeIdx idx, SPDirection d);
bool out_cor_ZX(SPNodeIdx idx, SPDirection d);

SPDirection commonDirection(SPNodeIdx idx1, SPDirection d);
SPNodeBoundary commonBoundary_XY(SPNodeIdx idx1, SPNodeIdx idx2);
SPDirection commonDirection_XY(SPNodeIdx idx1, SPNodeIdx idx2);
SPDirection commonDirection_XY(SPNodeIdx idx1, SPDirection d);

SPNodeBoundary commonBoundary_YZ(SPNodeIdx idx1, SPNodeIdx idx2);
SPDirection commonDirection_YZ(SPNodeIdx idx1, SPNodeIdx idx2);
SPDirection commonDirection_YZ(SPNodeIdx idx1, SPDirection d);

SPDirection commonDirection_ZX(SPNodeIdx idx1, SPDirection d);

SPDirection transfer_3Ddir_to_2Ddir(SPNodeIdx idx, SPDirection d);

int Direction_Decompose(const SPDirection&, SPDirection&, SPDirection&,
		SPDirection&);

SPDirection Direction_Compose(const SPDirection&, const SPDirection&); //2D

inline SPDirection oppositeDirection(const SPDirection& d) {
	return d == ErrSPDirection ? ErrSPDirection : SP_DIR_OP[int(d)];
}
inline bool isPlus(const SPDirection& d){
	ASSERT(d!=ErrSPDirection);
	return SP_DIR_MorP[int(d)]==1;
}
inline bool isMinus(const SPDirection& d){
	ASSERT(d!=ErrSPDirection);
	return SP_DIR_MorP[int(d)]==-1;
}
//  == to ==
SPNodeIdx toSPNodeIdx(int i);
SPNodeIdx toSPNodeIdx(SPDirection i);

SPNodeBoundary toSPNodeBoundary(int i);
SPNodeBoundary toSPNodeBoundary(SPDirection i);

SPDirection toDirection(int i);
SPDirection toDirection(SPNodeIdx i);
SPDirection toDirection(SPNodeBoundary i);

bool hasType(int type, SPNodeType t);
std::string parseQTNodeType(int type);
//end Function <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

template<typename Cell, typename Data, int Dim>
class SPNode: public ObjectBase {
public:
	typedef SPNode<Cell, Data, Dim> Node;
	typedef SPNode<Cell, Data, Dim>* pNode;
	typedef Cell* pCell;
	typedef Data* pData;
	typedef Cell Cell_type;
	typedef Data Data_type;
	typedef void (*pFun_process_SP)(pNode, utPointer);
	static const int DIM = Dim;
	static const int NUM_CELLS = Cell::NUM_VERTEXES;  //gcc 4.8
	static const int NUM_NEBOR = Cell::NUM_VERTEXES - 1;  //2D=3  3D=7
protected:
	int _nodetype;
	int _level;
	int _idx;

	inline int _height(const SPNode<Cell, Data, Dim>*) const;
	void _traversal(pNode pn, pFun_process_SP visit, utPointer p);
public:
	pCell cell;
	pNode father;
	pNode child[NUM_CELLS];
	pNode neighbor[NUM_NEBOR];

	pData data;

	SPNode(pNode f, int nt, int level, int idx, const Cell &c);
	~SPNode();

	inline int getType() const;
	inline void setType(int);

	inline int getLevel() const;
	inline int getIdx() const;

	inline bool hasChild() const;
	inline bool hasChild(int idx) const;
	inline bool isFullChild() const;
	inline int countChild() const;

	inline int Height() const;

	inline SPNodeIdx getEIdx() const;  //

	bool isAdjacent(SPNodeBoundary b) const;
	SPNodeIdx reflectIdx(SPDirection d) const; //Direction on Aix
	SPNodeIdx reflectIdx_Vertex() const;     //only for 3D
	bool isAdj_Vertex(SPDirection d) const;  //only for 3D

	inline void setNeighbors( //
			pNode, pNode, pNode, //
			pNode = NULL_PTR, pNode = NULL_PTR, //
			pNode = NULL_PTR, pNode = NULL_PTR);

	int creatChild_all(size_t max_level);
	void Traversal(pFun_process_SP, utPointer);
	int whichChild(const typename Cell::Point &p) const;  //revisable

	//shape relate function =====================
	typename Cell_type::Point getPoint(SPDirection d) const;
	//get =======================================

	//Data-process===============================
	bool isDataEmpty() const;
	//IO ========================================
	void show() const;
	void draw_to_vtu(std::string filename) const;
	void draw_to_gnuplot(std::string filename, int mode, SPPlane = SP_XY) const; //2D
	//Neighbor Finding============================
	pNode getAdjNeighborFast(SPDirection d);
	pNode getCorNeighborFast(SPDirection d);
	pNode getNeighbor_Adj_CorFast(SPDirection d);
	pNode getCorNeighbor_XYZFast(SPDirection d);
	pNode getNeighborFast(SPDirection d);
};

template<typename Cell, typename Data, int Dim>
SPNode<Cell, Data, Dim>::SPNode(pNode f, int nt, int level, int idx,
		const Cell &c) {
	_nodetype = nt;
	_level = level;
	cell = new Cell(c);
	father = f;
	_idx = idx;
	//_hasChild = false;
	data = NULL_PTR;
	for (int i = 0; i < this->NUM_CELLS; i++) {
		child[i] = NULL_PTR;
	}
	for (int i = 0; i < this->NUM_NEBOR; i++) {
		neighbor[i] = NULL_PTR;
	}
}
template<typename Cell, typename Data, int Dim>
SPNode<Cell, Data, Dim>::~SPNode() {
	delete cell;
	if (data != NULL_PTR) {
		delete data;
	}
}

template<typename Cell, typename Data, int Dim>
inline int SPNode<Cell, Data, Dim>::getType() const {
	return _nodetype;
}

template<typename Cell, typename Data, int Dim>
inline void SPNode<Cell, Data, Dim>::setType(int type) {
	_nodetype = type;
}

template<typename Cell, typename Data, int Dim>
inline int SPNode<Cell, Data, Dim>::getLevel() const {
	return _level;
}
template<typename Cell, typename Data, int Dim>
inline int SPNode<Cell, Data, Dim>::getIdx() const {
	return _idx;
}
template<typename Cell, typename Data, int Dim>
inline bool SPNode<Cell, Data, Dim>::hasChild() const {
	for (int i = 0; i < this->NUM_CELLS; i++) {
		if (this->child[i] != NULL_PTR) {
			return true;
		}
	}
	return false;
}
template<typename Cell, typename Data, int Dim>
inline bool SPNode<Cell, Data, Dim>::hasChild(int idx) const {
	return this->child[idx] != NULL_PTR;
}
template<typename Cell, typename Data, int Dim>
inline bool SPNode<Cell, Data, Dim>::isFullChild() const {
	bool res = this->child[0] != NULL_PTR;
	for (int i = 1; i < this->NUM_CELLS; i++) {
		res = res && (this->child[i] != NULL_PTR);
	}
	return res;
}
template<typename Cell, typename Data, int Dim>
inline int SPNode<Cell, Data, Dim>::countChild() const {
	int res = 0;
	for (int i = 0; i < this->NUM_CELLS; i++) {
		res += (this->child[i] != NULL_PTR) ? 1 : 0;
	}
	return res;
}

template<typename Cell, typename Data, int Dim>
inline SPNodeIdx SPNode<Cell, Data, Dim>::getEIdx() const {
	return SP_NODEIDX[_idx];
}
template<typename Cell, typename Data, int Dim>
int SPNode<Cell, Data, Dim>::_height(
		const SPNode<Cell, Data, Dim>* Current) const {
	if (Current == NULL_PTR) {
		return 0;
	}
	if (!Current->hasChild())
		return 0;
	else if (Dim == 2) {
		return 1
				+ MAX(_height(Current->child[0]), _height(Current->child[1]),
						_height(Current->child[2]), _height(Current->child[3]));
	} else { //Dim==3
		return 1
				+ MAX(_height(Current->child[0]), _height(Current->child[1]),
						_height(Current->child[2]), _height(Current->child[3]),
						_height(Current->child[4]), _height(Current->child[5]),
						_height(Current->child[6]), _height(Current->child[7]));
	}
}
template<typename Cell, typename Data, int Dim>
inline int SPNode<Cell, Data, Dim>::Height() const {
	return this->_height(this);
}

template<typename Cell, typename Data, int Dim>
bool SPNode<Cell, Data, Dim>::isAdjacent(SPNodeBoundary b) const {
	return b != ErrSPBoundary ? SP_ADJ[_idx][int(b) - 4] : false;
}
template<typename Cell, typename Data, int Dim>
SPNodeIdx SPNode<Cell, Data, Dim>::reflectIdx(SPDirection d) const {
	return (d >= 4 && d <= 9) ? toSPNodeIdx(SP_REFLECT[_idx][d - 4]) : ErrSPIdx;
}
template<typename Cell, typename Data, int Dim>
SPNodeIdx SPNode<Cell, Data, Dim>::reflectIdx_Vertex() const {
	return toSPNodeIdx(SP_DIANODE_XYZ[_idx]);
}
template<typename Cell, typename Data, int Dim>
bool SPNode<Cell, Data, Dim>::isAdj_Vertex(SPDirection d) const {
	// d 18+0 to 18+7
	if (d < 18 || d > 25) {
		return false;
	} else {
		return SP_VERTEX_XYZ[_idx][d - 18];
	}
}
template<typename Cell, typename Data, int Dim>
inline void SPNode<Cell, Data, Dim>::setNeighbors(
		//
		SPNode<Cell, Data, Dim>::pNode p1, SPNode<Cell, Data, Dim>::pNode p2,
		SPNode<Cell, Data, Dim>::pNode p3, SPNode<Cell, Data, Dim>::pNode p4,
		SPNode<Cell, Data, Dim>::pNode p5, SPNode<Cell, Data, Dim>::pNode p6,
		SPNode<Cell, Data, Dim>::pNode p7) {

	if (Dim == 3) {
		neighbor[0] = p1;
		neighbor[1] = p2;
		neighbor[2] = p3;
		neighbor[3] = p4;
		neighbor[4] = p5;
		neighbor[5] = p6;
		neighbor[6] = p7;
	} else {
		neighbor[0] = p1;
		neighbor[1] = p2;
		neighbor[2] = p4;
	}
}
template<typename Cell, typename Data, int Dim>
int SPNode<Cell, Data, Dim>::creatChild_all(size_t max_level) {
	size_t le = this->getLevel();
	if (this->hasChild() == false && le < max_level) {
		int ltmp = le + 1;
		this->child[0] = new Node(this, SPT_normal, ltmp, NWB,
				Cell(this->cell->getPoint(eCPL_M, eCPL_C, eCPL_M),
						this->cell->getPoint(eCPL_C, eCPL_P, eCPL_C)));
		this->child[0]->father = this;
		this->child[1] = new Node(this, SPT_normal, ltmp, SWB,
				Cell(this->cell->getPoint(eCPL_M, eCPL_M, eCPL_M),
						this->cell->getPoint(eCPL_C, eCPL_C, eCPL_C)));
		this->child[1]->father = this;
		this->child[2] = new Node(this, SPT_normal, ltmp, NEB,
				Cell(this->cell->getPoint(eCPL_C, eCPL_C, eCPL_M),
						this->cell->getPoint(eCPL_P, eCPL_P, eCPL_C)));
		this->child[2]->father = this;
		this->child[3] = new Node(this, SPT_normal, ltmp, SEB,
				Cell(this->cell->getPoint(eCPL_C, eCPL_M, eCPL_M),
						this->cell->getPoint(eCPL_P, eCPL_C, eCPL_C)));
		this->child[3]->father = this;
		if (Dim == 3) {
			this->child[4] = new Node(this, SPT_normal, ltmp, NWF,
					Cell(this->cell->getPoint(eCPL_M, eCPL_C, eCPL_C),
							this->cell->getPoint(eCPL_C, eCPL_P, eCPL_P)));
			this->child[4]->father = this;
			this->child[5] = new Node(this, SPT_normal, ltmp, SWF,
					Cell(this->cell->getPoint(eCPL_M, eCPL_M, eCPL_C),
							this->cell->getPoint(eCPL_C, eCPL_C, eCPL_P)));
			this->child[5]->father = this;
			this->child[6] = new Node(this, SPT_normal, ltmp, NEF,
					Cell(this->cell->getPoint(eCPL_C, eCPL_C, eCPL_C),
							this->cell->getPoint(eCPL_P, eCPL_P, eCPL_P)));
			this->child[6]->father = this;
			this->child[7] = new Node(this, SPT_normal, ltmp, SEF,
					Cell(this->cell->getPoint(eCPL_C, eCPL_M, eCPL_C),
							this->cell->getPoint(eCPL_P, eCPL_C, eCPL_P)));
			this->child[7]->father = this;
		}
		return 1;
	}
	return 0;
}
template<typename Cell, typename Data, int Dim>
void SPNode<Cell, Data, Dim>::_traversal(SPNode<Cell, Data, Dim>::pNode pn,
		SPNode<Cell, Data, Dim>::pFun_process_SP visit, utPointer p) {
	if (pn == NULL_PTR) {
		return;
	} else {
		(*visit)(pn, p);
		if (pn->hasChild()) {
			for (int i = 0; i < NUM_CELLS; i++) {
				pNode c = pn->child[i];
				if (c != NULL_PTR) {
					_traversal(c, visit, p);
				}
			}
		}
	}
}

template<typename Cell, typename Data, int Dim>
void SPNode<Cell, Data, Dim>::Traversal(
		SPNode<Cell, Data, Dim>::pFun_process_SP visit, utPointer p) {
	_traversal(this, visit, p);
}
template<typename Cell, typename Data, int Dim>
int SPNode<Cell, Data, Dim>::whichChild(const typename Cell::Point &p) const {
	typename Cell::Point mm = cell->getPoint(eCPL_M, eCPL_M, eCPL_M);
	typename Cell::Point pp = cell->getPoint(eCPL_P, eCPL_P, eCPL_P);
	if (p.x < mm.x || p.x > pp.x || p.y < mm.y || p.y > pp.y
			|| (Dim == 3 ? (p[2] < mm[2] || p[2] > pp[2]) : false)) {
		return -1;
	}
	mm = cell->getPoint(eCPL_M, eCPL_C, eCPL_M);
	pp = cell->getPoint(eCPL_C, eCPL_P, eCPL_C);
	if (p.x >= mm.x && p.x <= pp.x && p.y >= mm.y && p.y <= pp.y
			&& (Dim == 3 ? (p[2] >= mm[2] && p[2] <= pp[2]) : true)) {
		return 0;
	}
	mm = cell->getPoint(eCPL_M, eCPL_M, eCPL_M);
	pp = cell->getPoint(eCPL_C, eCPL_C, eCPL_C);
	if (p.x >= mm.x && p.x <= pp.x && p.y >= mm.y && p.y <= pp.y
			&& (Dim == 3 ? (p[2] >= mm[2] && p[2] <= pp[2]) : true)) {
		return 1;
	}
	mm = cell->getPoint(eCPL_C, eCPL_C, eCPL_M);
	pp = cell->getPoint(eCPL_P, eCPL_P, eCPL_C);
	if (p.x >= mm.x && p.x <= pp.x && p.y >= mm.y && p.y <= pp.y
			&& (Dim == 3 ? (p[2] >= mm[2] && p[2] <= pp[2]) : true)) {
		return 2;
	}
	mm = cell->getPoint(eCPL_C, eCPL_M, eCPL_M);
	pp = cell->getPoint(eCPL_P, eCPL_C, eCPL_C);
	if (p.x >= mm.x && p.x <= pp.x && p.y >= mm.y && p.y <= pp.y
			&& (Dim == 3 ? (p[2] >= mm[2] && p[2] <= pp[2]) : true)) {
		return 3;
	}
	if (Dim == 3) {
		mm = cell->getPoint(eCPL_M, eCPL_C, eCPL_C);
		pp = cell->getPoint(eCPL_C, eCPL_P, eCPL_P);
		if (p.x >= mm.x && p.x <= pp.x && p.y >= mm.y && p.y <= pp.y
				&& p[2] >= mm[2] && p[2] <= pp[2]) {
			return 4;
		}
		mm = cell->getPoint(eCPL_M, eCPL_M, eCPL_C);
		pp = cell->getPoint(eCPL_C, eCPL_C, eCPL_P);
		if (p.x >= mm.x && p.x <= pp.x && p.y >= mm.y && p.y <= pp.y
				&& p[2] >= mm[2] && p[2] <= pp[2]) {
			return 5;
		}
		mm = cell->getPoint(eCPL_C, eCPL_C, eCPL_C);
		pp = cell->getPoint(eCPL_P, eCPL_P, eCPL_P);
		if (p.x >= mm.x && p.x <= pp.x && p.y >= mm.y && p.y <= pp.y
				&& p[2] >= mm[2] && p[2] <= pp[2]) {
			return 6;
		}
		mm = cell->getPoint(eCPL_C, eCPL_M, eCPL_C);
		pp = cell->getPoint(eCPL_P, eCPL_C, eCPL_P);
		if (p.x >= mm.x && p.x <= pp.x && p.y >= mm.y && p.y <= pp.y
				&& p[2] >= mm[2] && p[2] <= pp[2]) {
			return 7;
		}
	}
	return -1;
}
const CellPointLocation Dir_to_CPL[26][3] = { //
											  //         X         y      z
		{ eCPL_M, eCPL_P, eCPL_C }, //SPD_MP = 0,
				{ eCPL_M, eCPL_M, eCPL_C }, //SPD_MM = 1,
				{ eCPL_P, eCPL_P, eCPL_C }, //SPD_PP = 2,
				{ eCPL_P, eCPL_M, eCPL_C }, //SPD_PM = 3,
				{ eCPL_M, eCPL_C, eCPL_C },	//SPD_IM = 4,
				{ eCPL_C, eCPL_P, eCPL_C },	//SPD_JP = 5,
				{ eCPL_P, eCPL_C, eCPL_C },	//SPD_IP = 6,
				{ eCPL_C, eCPL_M, eCPL_C },	//SPD_JM = 7,
				{ eCPL_C, eCPL_C, eCPL_P },	//SPD_KP = 8,
				{ eCPL_C, eCPL_C, eCPL_M },	//SPD_KM = 9,
				{ eCPL_C, eCPL_M, eCPL_P }, //SPD_MP = 10 + 0,
				{ eCPL_C, eCPL_M, eCPL_M }, //SPD_MM = 10 + 1
				{ eCPL_C, eCPL_P, eCPL_P }, //SPD_PP = 10 + 2,
				{ eCPL_C, eCPL_P, eCPL_M }, //SPD_PM = 10 + 3,
				{ eCPL_P, eCPL_C, eCPL_M }, //SPD_MP = 14 + 0,
				{ eCPL_M, eCPL_C, eCPL_M }, //SPD_MM = 14 + 1,
				{ eCPL_P, eCPL_C, eCPL_P }, //SPD_PP = 14 + 2,
				{ eCPL_M, eCPL_C, eCPL_P }, //SPD_PM = 14 + 3,
				{ eCPL_M, eCPL_M, eCPL_M }, //SPD_MMM = 17 + 1,
				{ eCPL_P, eCPL_M, eCPL_M }, //SPD_PMM = 17 + 2,
				{ eCPL_M, eCPL_P, eCPL_M }, //SPD_MPM = 17 + 3,
				{ eCPL_P, eCPL_P, eCPL_M }, //SPD_PPM = 17 + 4,
				{ eCPL_M, eCPL_M, eCPL_P }, //SPD_MMP = 17 + 5,
				{ eCPL_P, eCPL_M, eCPL_P }, //SPD_PMP = 17 + 6,
				{ eCPL_M, eCPL_P, eCPL_P }, //SPD_MPP = 17 + 7,
				{ eCPL_P, eCPL_P, eCPL_P }  //SPD_PPP = 17 + 8,
		};

const SPDirection FaceDir_to_VertexDir2D[4][2] = { //
												   //         X         y
		{ SPD_MM, SPD_MP },	//SPD_IM = 4,
				{ SPD_MP, SPD_PP },	//SPD_JP = 5,
				{ SPD_PM, SPD_PP },	//SPD_IP = 6,
				{ SPD_MM, SPD_PM },	//SPD_JM = 7,
		};

template<typename Cell, typename Data, int Dim>
typename SPNode<Cell, Data, Dim>::Cell_type::Point SPNode<Cell, Data, Dim>::getPoint(
		SPDirection d) const {
	ASSERT(d != ErrSPDirection);
	if (Dim == 2) {
		ASSERT(d <= 7);
		return this->cell->getPoint(Dir_to_CPL[int(d)][0],
				Dir_to_CPL[int(d)][1], Dir_to_CPL[int(d)][2]);
	} else { //Dim == 3
		return this->cell->getPoint(Dir_to_CPL[int(d)][0],
				Dir_to_CPL[int(d)][1], Dir_to_CPL[int(d)][2]);
	}
}

template<typename Cell, typename Data, int Dim>
bool SPNode<Cell, Data, Dim>::isDataEmpty() const {
	if (data == NULL_PTR) {
		return true;
	} else {
		return false;
	}
}
template<typename Cell, typename Data, int Dim>
void SPNode<Cell, Data, Dim>::show() const {
	std::cout << "OCNode ------\n";
	if (NULL_PTR == this) {
		std::cout << "Empty Node \n";
		return;
	}
	std::cout << "Dimension = " << DIM << "D\n";
	std::cout << "level      :" << this->_level << "\n";
	std::cout << "node type  :" << this->_nodetype << "\n";
	std::cout << "idx        :" << this->_idx << "\n";
	std::cout << "CELL  show =========\n";
	std::cout << std::scientific;
	std::cout << "Dx         :" << this->cell->getDx() << "\n";
	std::cout << "Dx         :" << this->cell->getDy() << "\n";
	if (Dim == 3) {
		std::cout << "Dx         :" << this->cell->getDz() << "\n";
	}
	this->cell->show();
	std::cout << "DATA  show =========\n";
	if (NULL_PTR == this->data) {
		std::cout << "No data \n";
	} else {
		this->data->show_info();
	}
}
template<typename Cell, typename Data, int Dim>
void SPNode<Cell, Data, Dim>::draw_to_vtu(std::string filename) const {
	LarusDef::size_type dim = Dim;
	LarusDef::size_type vertexes = Cell::NUM_VERTEXES;
	std::ofstream fs;
	open_file(fs, filename, 1);
	int numnode = 1;
	vtu_unstructured_grid_file_head(fs);
	fs << "<Piece NumberOfPoints=\"" << numnode * vertexes
			<< "\" NumberOfCells=\"" << numnode << "\">";
	fs << "<Points>";
	fs
			<< "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";
	this->cell->output_vertex_in_vtk_order(fs);
	fs << "</DataArray>";
	fs << "</Points>";
	fs << "<Cells>";
	fs << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
	for (int i = 0; i < numnode * vertexes; i++) {
		fs << i << " ";
	}
	fs << "</DataArray>";
	fs << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
	for (int i = 1; i < numnode + 1; i++) {
		fs << i * vertexes << " ";
	}
	fs << "</DataArray>";
	fs << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
	for (int i = 0; i < numnode; i++) {
		if (dim == 3) {
			fs << 11 << " ";
		} else if (dim == 2) {
			fs << 8 << " ";
		}
	}
	fs << "</DataArray>";
	fs << "</Cells>";
	fs << "</Piece>";
	vtu_unstructured_grid_file_end(fs);
	fs.close();
}

//SPD_MP = 0,
//SPD_MM = 1,
//SPD_PP = 2,
//SPD_PM = 3,
const int NEI_ONE_STEP_ADJ_IDX[8][6] = {
//       4      5      6      7      8      9
//       SPD_IM SPD_JP SPD_IP SPD_JM SPD_KP SPD_KM
		{ 0, 1, -1, -1, -1, 2 }, //0 MPM
		{ 0, -1, -1, 1, -1, 2 }, //1 MMM
		{ -1, 1, 0, -1, -1, 2 }, //2 PPM
		{ -1, -1, 0, 1, -1, 2 }, //3 PMM
		{ 0, 1, -1, -1, 2, -1 }, //4 MPP
		{ 0, -1, -1, 1, 2, -1 }, //5 MMP
		{ -1, 1, 0, -1, 2, -1 }, //6 PPP
		{ -1, -1, 0, 1, 2, -1 }, //7 PMP
		};
template<typename Cell, typename Data, int Dim>
typename SPNode<Cell, Data, Dim>::pNode SPNode<Cell, Data, Dim>::getAdjNeighborFast(
		SPDirection d) {
	//d form 4 to 9
	//this is not empty
	//this has father
	if (NEI_ONE_STEP_ADJ_IDX[_idx][d - 4] != -1) {
		return neighbor[NEI_ONE_STEP_ADJ_IDX[_idx][d - 4]];
	} else {
		return father->child[this->reflectIdx(d)];
	}
}
const int IDX_YZ_TO_IDX_XY[8] = { 5, 1, 4, 0, 7, 3, 6, 2 };
const int IDX_ZX_TO_IDX_XY[8] = { 3, 1, 7, 5, 2, 0, 6, 4 };
//1  one step
//2  two step in
//-2 two step out
const int NEI_ONE_STEP_COR_IDX_XY[8][4] = {
//
		{ 1, -2, -2, 2 }, //0 MPM
		{ -2, 1, 2, -2 }, //1 MMM
		{ -2, 2, 1, -2 }, //2 PPM
		{ 2, -2, -2, 1 }, //3 PMM
		{ 1, -2, -2, 2 }, //4 MPP
		{ -2, 1, 2, -2 }, //5 MMP
		{ -2, 2, 1, -2 }, //6 PPP
		{ 2, -2, -2, 1 }, //7 PMP
		};
inline int cor_step(SPNodeIdx idx, SPDirection d) {
	if (d == ErrSPDirection)
		return -3;
	if (d >= 0 && d <= 3)
		return NEI_ONE_STEP_COR_IDX_XY[int(idx)][int(d)];
	if (d >= 10 && d <= 13)
		return NEI_ONE_STEP_COR_IDX_XY[IDX_YZ_TO_IDX_XY[int(idx)]][int(d) - 10];
	if (d >= 14 && d <= 17)
		return NEI_ONE_STEP_COR_IDX_XY[IDX_ZX_TO_IDX_XY[int(idx)]][int(d) - 14];
	return -3;
}

template<typename Cell, typename Data, int Dim>
typename SPNode<Cell, Data, Dim>::pNode SPNode<Cell, Data, Dim>::getCorNeighborFast(
		SPDirection d) {
	//this is not empty
	//this has father
	pNode ca = NULL_PTR;    //common ancestor
	if (father != NULL_PTR && cor_step(getEIdx(), d) != 2) {
		if (cor_step(getEIdx(), d) == 1) {
			if (Dim == 3) {
				if (d >= 0 && d <= 3)  //XY
					return ca = neighbor[3];
				if (d >= 10 && d <= 13)  //YZ
					return ca = neighbor[4];
				if (d >= 14 && d <= 17)  //ZX
					return ca = neighbor[5];
			} else { //Dim==2
				return ca = neighbor[2];
			}
		} else {
			ca = father->getAdjNeighborFast(commonDirection(getEIdx(), d));
		}
	} else {
		ca = father;
	}
	//Follow opposite path to locate the neighbor
	pNode pt = NULL_PTR;
	if (ca != NULL_PTR && ca->hasChild()) {
		pt = ca->child[diagonalIdx(getEIdx(), d)];
	} else {
		pt = ca;
	}
	return pt;
}
template<typename Cell, typename Data, int Dim>
typename SPNode<Cell, Data, Dim>::pNode SPNode<Cell, Data, Dim>::getNeighbor_Adj_CorFast(
		SPDirection d) {
	// AdjNeighbor
	if (d >= 4 && d <= 9) {
		return getAdjNeighborFast(d);
	}
	if (d != ErrSPDirection && d <= 18) {
		return getCorNeighborFast(d);
	}
	return NULL_PTR;
}
template<typename Cell, typename Data, int Dim>
typename SPNode<Cell, Data, Dim>::pNode SPNode<Cell, Data, Dim>::getCorNeighbor_XYZFast(
		SPDirection d) {
	// d  17+1 to 17+8
	pNode ca = NULL_PTR;    //common ancestor
	if (father != NULL_PTR && isAdj_Vertex(d)) {
		int Idx_xyz = 7;
		ca = neighbor[Idx_xyz];
	} else {
		SPDirection newd = transfer_3Ddir_to_2Ddir(getEIdx(), d);
		if (newd != ErrSPDirection) {
			ca = father->getNeighbor_Adj_CorFast(newd);
		} else {
			ca = father;
		}
	}
	pNode pt = NULL_PTR;
	if (ca != NULL_PTR && ca->hasChild()) {
		pt = ca->child[reflectIdx_Vertex()];
	} else {
		pt = ca;
	}
	return pt;
}

template<typename Cell, typename Data, int Dim>
typename SPNode<Cell, Data, Dim>::pNode SPNode<Cell, Data, Dim>::getNeighborFast(
		SPDirection d) {
// AdjNeighbor
	if (d != ErrSPDirection && d <= 18) {
		return getNeighbor_Adj_CorFast(d);
	} else if (d >= 18 && d <= 25) {
		return getCorNeighbor_XYZFast(d);
	}
	return NULL_PTR;
}

//===============================================
//end of SPNode==================================
//===============================================

//QTNode=========================================
typedef SPNode<Cell2D, CellData2D, 2> QTNode;
typedef QTNode* pQTNode;

//OCNode=========================================
typedef SPNode<Cell3D, CellData3D, 3> OCNode;
typedef OCNode* pOCNode;

//Function out of class==========================
//================================================
pQTNode getChild_vertex(pQTNode pn, SPDirection dir);

template<typename NODE>
NODE* getFirstChild(const NODE* p) {
	NODE* c = NULL_PTR;
	for (int i = 0; i < NODE::NUM_CELLS; ++i) {
		c = p->child[i];
		if (c != NULL_PTR) {
			return c;
		}
	}
	return c;
}

template<typename NODE>
NODE* getLastChild(const NODE* p) {
	NODE* c = NULL_PTR;
	for (int i = NODE::NUM_CELLS - 1; i > 0; --i) {
		c = p->child[i];
		if (c != NULL_PTR) {
			return c;
		}
	}
	return c;
}

template<typename NODE>
NODE* getSiblingPlus(NODE* p) {
	NODE* f = p->father;
	if (f == NULL_PTR) {
		return p;
	}
	for (int i = p->getIdx() + 1; i < NODE::NUM_CELLS; ++i) {
		NODE* c = f->child[i];
		if (c != NULL_PTR) {
			return c;
		}
	}
	return getSiblingPlus(f);
}

template<typename NODE>
NODE* getSiblingMinus(NODE* p) {
	NODE* f = p->father;
	if (f == NULL_PTR) {
		return p;
	}
	for (int i = p->getIdx() - 1; i >= 0; --i) {
		NODE* c = f->child[i];
		if (c != NULL_PTR) {
			return c;
		}
	}
	return getSiblingMinus(f);
}

template<typename NODE>
NODE* getFirstLeaf(NODE* p) {
	NODE* c = getFirstChild(p);
	if (c == NULL_PTR) {
		return p;
	} else {
		NODE* resc = c;
		while (c != NULL_PTR) {
			resc = c;
			c = getFirstChild(c);
		}
		return resc;
	}
}

template<typename NODE>
NODE* getLastLeaf(NODE* p) {
	NODE* c = getLastChild(p);
	if (c == NULL_PTR) {
		return p;
	} else {
		NODE* resc = c;
		while (c != NULL_PTR) {
			resc = c;
			c = getLastChild(c);
		}
		return resc;
	}
}

int point_at_which_child(pQTNode pnode, Point2D& p);
int point_at_which_child(pOCNode pnode, Point3D& p);

int relation_node_level(pQTNode, pQTNode);
int relation_node_level(pOCNode, pOCNode);

void draw_boundary(pQTNode qnode, FILE* file);

template<typename Cell, typename Data, int Dim>
void draw_gnuplot_boundary(const SPNode<Cell, Data, Dim>* pnode, FILE* file,
		SPPlane spp) {
	if (spp == ErrSPPlane) {
		return;
	}
	arrayListT<typename Cell::Point> listp(4);
	for (int i = 0; i < 4; i++) {
		listp[i] = pnode->getPoint(toDirection(int(spp) + i));
	}
	int x = (spp == SP_XY) ? 0 : (spp == SP_YZ) ? 1 : 2;
	int y = (spp == SP_XY) ? 1 : (spp == SP_YZ) ? 2 : 0;
	fprintf(file, "%f %f\n", listp[0][x], listp[0][y]);
	fprintf(file, "%f %f\n\n", listp[1][x], listp[1][y]);
	fprintf(file, "%f %f\n", listp[1][x], listp[1][y]);
	fprintf(file, "%f %f\n\n", listp[3][x], listp[3][y]);
	fprintf(file, "%f %f\n", listp[3][x], listp[3][y]);
	fprintf(file, "%f %f\n\n", listp[2][x], listp[2][y]);
	fprintf(file, "%f %f\n", listp[2][x], listp[2][y]);
	fprintf(file, "%f %f\n\n", listp[0][x], listp[0][y]);
}

template<typename Cell, typename Data, int Dim>
void SPNode<Cell, Data, Dim>::draw_to_gnuplot(std::string filename, int mode,
		SPPlane spp) const {    //2D
	FILE* data = open_file(filename, mode);
	draw_gnuplot_boundary(this, data, spp);
	fclose(data);
}

//===============================================
//===============================================

enum SPNodeFaceType {
	SPFT_Error = -1,
	SPFT_Boundary = 0,
	SPFT_Equal = 1,
	SPFT_FineCoarse = 2,
	SPFT_CoarseFine = 3,
};

std::string parseSPNodeFaceType(const SPNodeFaceType& t);

template<typename Cell, typename Data, int Dim>
SPNodeFaceType getFaceType(  //
		SPNode<Cell, Data, Dim>* p,     //main node
		SPNode<Cell, Data, Dim>* pn) {
	if (p == NULL_PTR) {
		return SPFT_Error;
	}
	if (pn == NULL_PTR) {
		return SPFT_Boundary;
	}
	if (pn->getType() == SPT_ghost) {
		return SPFT_Boundary;
	}
	if (p->getLevel() > pn->getLevel())
		return SPFT_FineCoarse;
	if (p->getLevel() == pn->getLevel() && !pn->hasChild())
		return SPFT_Equal;
	return SPFT_CoarseFine;
}

template<typename Cell, typename Data, int Dim>
class SPNodeFace: public ObjectBase {
public:
	typedef SPNode<Cell, Data, Dim> Node;
	typedef SPNode<Cell, Data, Dim>* pNode;
public:
	pNode pnode;
	pNode pneighbor;
	SPNodeFaceType face_type;
	SPDirection direction;

	SPNodeFace();
	SPNodeFace(pNode, pNode, SPDirection, SPNodeFaceType);

	SPNodeFace(const SPNodeFace&);
	//
	//operator ==================================
	SPNodeFace& operator=(const SPNodeFace &a);
	bool operator==(const SPNodeFace& a) const {
		return (pnode == a.pnode) && (pneighbor == a.pneighbor)
				&& (face_type == a.face_type) && (direction == a.direction);
	}
	bool operator!=(const SPNodeFace& a) const {
		return !((pnode == a.pnode) && (pneighbor == a.pneighbor)
				&& (face_type == a.face_type) && (direction == a.direction));
	}
	//show=======================================
	void show() const;
};
template<typename Cell, typename Data, int Dim>
SPNodeFace<Cell, Data, Dim>::SPNodeFace() {
	pnode = NULL_PTR;
	pneighbor = NULL_PTR;
	face_type = SPFT_Error;
	direction = ErrSPDirection;
}

template<typename Cell, typename Data, int Dim>
SPNodeFace<Cell, Data, Dim>::SPNodeFace(
		typename SPNodeFace<Cell, Data, Dim>::pNode p,  //
		typename SPNodeFace<Cell, Data, Dim>::pNode pn,  //
		SPDirection dir,  //
		SPNodeFaceType ft) { //
	pnode = p;
	pneighbor = pn;
	face_type = ft;
	direction = dir;
}
template<typename Cell, typename Data, int Dim>
SPNodeFace<Cell, Data, Dim>::SPNodeFace(const SPNodeFace<Cell, Data, Dim>& a) {
	pnode = a.pnode;
	pneighbor = a.pneighbor;
	face_type = a.face_type;
	direction = a.direction;
}
template<typename Cell, typename Data, int Dim>
SPNodeFace<Cell, Data, Dim>& SPNodeFace<Cell, Data, Dim>::operator=(
		const SPNodeFace<Cell, Data, Dim> &a) {
	pnode = a.pnode;
	pneighbor = a.pneighbor;
	face_type = a.face_type;
	direction = a.direction;
	return *this;
}
template<typename Cell, typename Data, int Dim>
void SPNodeFace<Cell, Data, Dim>::show() const {
	std::cout << "Face show ===============\n";
	std::cout << "Dim        : " << Dim << '\n';
	std::cout << "Direction  : " << direction << '\n';
	std::cout << "Face type  : " << parseSPNodeFaceType(face_type) << '\n';
	std::cout << "Center p x : " << pnode->cell->getX(eCPL_C) << '\n';
	std::cout << "       p y : " << pnode->cell->getY(eCPL_C) << '\n';
	if (Dim == 3) {
		std::cout << "       p z : " << pnode->cell->getZ(eCPL_C) << '\n';
	}
	if (Dim == 2) {

	} else {

	}
}
//QTNodeFace=========================================
typedef SPNodeFace<Cell2D, CellData2D, 2> QTNodeFace;
typedef QTNodeFace* pQTNodeFace;

//OCNodeFace=========================================
typedef SPNodeFace<Cell3D, CellData3D, 3> OCNodeFace;
typedef OCNodeFace* pOCNodeFace;

//---------------------------------------------------
inline pQTNode new_ghost_node(const QTNodeFace& face) {
	int ghost_node_type = SPT_ghost;
	//direction 4 5 6 7
	Cell2D c(*(face.pnode->cell));
	switch (face.direction) {
	case SPD_IM:
		c.transfer(-c.getDx(), 0);
		break;
	case SPD_JM:
		c.transfer(0, -c.getDy());
		break;
	case SPD_IP:
		c.transfer(c.getDx(), 0);
		break;
	case SPD_JP:
		c.transfer(0, c.getDx());
		break;
	default:
		ASSERT(false);
	}
	pQTNode ghostnode = new QTNode(NULL_PTR, ghost_node_type,
			face.pnode->getLevel(), face.pnode->getIdx(), c);
	ghostnode->father = face.pnode;
	ghostnode->data = new QTNode::Data_type((*face.pnode->data));
	return ghostnode;
}

inline void delete_ghost_node(QTNodeFace& face) {
	if (face.pneighbor == NULL_PTR || face.pneighbor->getType() != SPT_ghost) {
		return;
	}
	delete face.pneighbor;
	face.pneighbor = NULL_PTR;
}

//
template<typename Cell, typename Data, int Dim>
class SPNodeVertex: public ObjectBase {
public:
	typedef SPNode<Cell, Data, Dim> Node;
	typedef SPNode<Cell, Data, Dim>* pNode;
public:
	pNode pnode;
	SPDirection direction;

	SPNodeVertex();
	SPNodeVertex(pNode, SPDirection);

	SPNodeVertex(const SPNodeVertex&);
//
//operator ==================================
	SPNodeVertex& operator=(const SPNodeVertex &a);
	bool operator==(const SPNodeVertex& a) const {
		return (pnode == a.pnode) && (direction == a.direction);
	}
	bool operator!=(const SPNodeVertex& a) const {
		return !((pnode == a.pnode) && (direction == a.direction));
	}
};
template<typename Cell, typename Data, int Dim>
SPNodeVertex<Cell, Data, Dim>::SPNodeVertex() {
	pnode = NULL_PTR;
	direction = ErrSPDirection;
}
template<typename Cell, typename Data, int Dim>
SPNodeVertex<Cell, Data, Dim>::SPNodeVertex(pNode pn, SPDirection dir) {
	pnode = pn;
	direction = dir;
}
template<typename Cell, typename Data, int Dim>
SPNodeVertex<Cell, Data, Dim>::SPNodeVertex(const SPNodeVertex& a) {
	pnode = a.pnode;
	direction = a.direction;
}
template<typename Cell, typename Data, int Dim>
SPNodeVertex<Cell, Data, Dim>& SPNodeVertex<Cell, Data, Dim>::operator=(
		const SPNodeVertex<Cell, Data, Dim> &a) {
	pnode = a.pnode;
	direction = a.direction;
	return *this;
}
//QTNodeVertex=========================================
typedef SPNodeVertex<Cell2D, CellData2D, 2> QTNodeVertex;
typedef QTNodeVertex* pQTNodeVertex;

//OCNodeVertex=========================================
typedef SPNodeVertex<Cell3D, CellData3D, 3> OCNodeVertex;
typedef OCNodeVertex* pOCNodeVertex;

//SPNode_iterator====================================
template<class Cell, class Data, int Dim, class _Ref, class _Ptr>
class _SPNode_iterator {
public:
	typedef LarusDef::size_type size_type;
	typedef LarusDef::size_type difference_type;
	typedef bidirectional_iterator_tag iterator_category;

	typedef SPNode<Cell, Data, Dim> _Node;

	typedef _SPNode_iterator<Cell, Data, Dim, _Node&, _Node*> iterator;
	typedef _SPNode_iterator<Cell, Data, Dim, const _Node&, const _Node*> const_iterator;
	typedef _SPNode_iterator<Cell, Data, Dim, _Ref, _Ptr> _Self;

	typedef _Node value_type;
	typedef _Ptr pointer;
	typedef _Ref reference;

	_Node* _ptr;

	_SPNode_iterator() {
		_ptr = NULL_PTR;
	}
	_SPNode_iterator(_Node* _x) {
		this->_ptr = _x;
	}
	_SPNode_iterator(const iterator& _x) {
		this->_ptr = _x._ptr;
	}

	void _incr() {
		if (_ptr->father != NULL_PTR) {
			_Node* s = getSiblingPlus(_ptr);
			if (s->father != NULL_PTR) {
				_ptr = getFirstLeaf(s);
			} else {
				_ptr = s;
			}
		} else {
			_ptr = getFirstLeaf(_ptr);
		}
	}

	void _decr() {
		if (_ptr->father != NULL_PTR) {
			_Node* s = getSiblingMinus(_ptr);
			if (s->father != NULL_PTR) {
				_ptr = getLastLeaf(s);
			} else {
				_ptr = s;
			}
		} else {
			_ptr = getLastLeaf(_ptr);
		}
	}

	bool operator==(const _SPNode_iterator& _x) const {
		return _ptr == _x._ptr;
	}
	bool operator!=(const _SPNode_iterator& _x) const {
		return _ptr != _x._ptr;
	}

	reference operator*() const {
		return (*_ptr);
	}

	pointer operator->() const {
		return &(operator*());
	}

	_Self & operator++() {
		this->_incr();
		return *this;
	}

	_Self operator++(int) {
		_Self __tmp = *this;
		this->_incr();
		return __tmp;
	}

	_Self & operator--() {
		this->_decr();
		return *this;
	}

	_Self operator--(int) {
		_Self __tmp = *this;
		this->_decr();
		return __tmp;
	}

	bool isExist() {
		return _ptr != NULL_PTR;
	}
};

} //this is the end of namespace

#endif /* QUADTREENODE_H_ */
