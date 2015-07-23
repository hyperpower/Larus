/************************
 //  \file   CellData.cpp
 //  \brief
 // 
 //  \author zhou
 //  \date   13 sept. 2014 
 ***********************/

#include "CellData.h"
#include <iostream>

namespace Larus {

//Class CellData2D ==================

CellData2D::CellData2D() :
		aCenterData() {
	utp_data = NULL_PTR;
	//aFaceData[0] = Array();
	//aFaceData[1] = Array();
	//aFaceData[2] = Array();
	//aFaceData[3] = Array();
}

CellData2D::CellData2D(int lc, int l1, int l2, int l3, int l4) :
		aCenterData(lc) {
	utp_data = NULL_PTR;
	aFaceData[0] = arrayList(l1);
	aFaceData[1] = arrayList(l2);
	aFaceData[2] = arrayList(l3);
	aFaceData[3] = arrayList(l4);
}

CellData2D::CellData2D(CellData2D &a) {
	utp_data = NULL_PTR;
	aCenterData = a.aCenterData;
	aFaceData[0] = a.aFaceData[0];
	aFaceData[1] = a.aFaceData[1];
	aFaceData[2] = a.aFaceData[2];
	aFaceData[3] = a.aFaceData[3];
}

CellData2D& CellData2D::operator=(const CellData2D &a) {
	if (this == &a) {
		return *this;
	} else {
		aCenterData = a.aCenterData;
		aFaceData[0] = a.aFaceData[0];
		aFaceData[1] = a.aFaceData[1];
		aFaceData[2] = a.aFaceData[2];
		aFaceData[3] = a.aFaceData[3];
		return *this;
	}
}

void CellData2D::reconstruct(int lc, int l1, int l2, int l3, int l4) {
	aCenterData.reconstruct(lc);
	aFaceData[0].reconstruct(l1);
	aFaceData[1].reconstruct(l2);
	aFaceData[2].reconstruct(l3);
	aFaceData[3].reconstruct(l4);
}

bool CellData2D::isEmpty() const {
	if (aCenterData.size() == 0 && aFaceData[0].size() == 0
			&& aFaceData[1].size() == 0 && aFaceData[2].size() == 0
			&& aFaceData[3].size() == 0) {
		return true;
	} else {
		return false;
	}
}

CellData2D::~CellData2D() {
}

void CellData2D::show_info() const {
	if (NULL_PTR == this) {
		std::cout << "No data \n";
	} else {
		std::cout << "center data:" << this->aCenterData.size() << std::endl;
		std::cout << "face data   :" << std::endl;
		std::cout << "  0         :" << this->aFaceData[0].size() << std::endl;
		std::cout << "  1         :" << this->aFaceData[1].size() << std::endl;
		std::cout << "  2         :" << this->aFaceData[2].size() << std::endl;
		std::cout << "  3         :" << this->aFaceData[3].size() << std::endl;
		if (!this->aCenterData.empty()) {
			std::cout << "CenterData=============>>" << std::endl;
			this->aCenterData.show();
			std::cout << "<<=======================" << std::endl;
		}
	}
}

//=================================

//Class CellData3D=================
CellData3D::CellData3D() :
		aCenterData() {
	utp_data = NULL_PTR;
	//aFaceData[0] = Array();
	//aFaceData[1] = Array();
	//aFaceData[2] = Array();
	//aFaceData[3] = Array();
	//aFaceData[4] = Array();
	//aFaceData[5] = Array();
}

CellData3D::CellData3D(int lc, int l1, int l2, int l3, int l4, int l5, int l6) :
		aCenterData(lc) {
	utp_data = NULL_PTR;
	aFaceData[0] = arrayList(l1);
	aFaceData[1] = arrayList(l2);
	aFaceData[2] = arrayList(l3);
	aFaceData[3] = arrayList(l4);
}

CellData3D::CellData3D(CellData3D &a) {
	utp_data = NULL_PTR;
	aCenterData = a.aCenterData;
	aFaceData[0] = a.aFaceData[0];
	aFaceData[1] = a.aFaceData[1];
	aFaceData[2] = a.aFaceData[2];
	aFaceData[3] = a.aFaceData[3];
	aFaceData[4] = a.aFaceData[4];
	aFaceData[5] = a.aFaceData[5];

}

CellData3D& CellData3D::operator=(const CellData3D &a) {
	if (this == &a) {
		return *this;
	} else {
		aCenterData = a.aCenterData;
		aFaceData[0] = a.aFaceData[0];
		aFaceData[1] = a.aFaceData[1];
		aFaceData[2] = a.aFaceData[2];
		aFaceData[3] = a.aFaceData[3];
		aFaceData[4] = a.aFaceData[4];
		aFaceData[5] = a.aFaceData[5];
		return *this;
	}
}

CellData3D::~CellData3D() {

}

void CellData3D::show_info() const {
	if (NULL_PTR == this) {
		std::cout << "No data \n";
	} else {
		std::cout << "center data:" << this->aCenterData.size() << std::endl;
		std::cout << "face data   :" << std::endl;
		std::cout << "  0         :" << this->aFaceData[0].size() << std::endl;
		std::cout << "  1         :" << this->aFaceData[1].size() << std::endl;
		std::cout << "  2         :" << this->aFaceData[2].size() << std::endl;
		std::cout << "  3         :" << this->aFaceData[3].size() << std::endl;
		std::cout << "  4         :" << this->aFaceData[4].size() << std::endl;
		std::cout << "  5         :" << this->aFaceData[5].size() << std::endl;
		if (!this->aCenterData.empty()) {
			std::cout << "=======================>>" << std::endl;
			this->aCenterData.show();
			std::cout << "<<=======================" << std::endl;
		}
	}

}

} //this is end of namespace
