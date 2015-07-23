/************************
 //  \file   Celldata.h
 //  \brief
 // 
 //  \author zhou
 //  \date   13 sept. 2014 
 ***********************/
#ifndef CELLDATA_H_
#define CELLDATA_H_

#include "../TypeDef.h"
#include "../Utility/ArrayList.h"


namespace Larus{

    class CellData2D : public ObjectBase{

    public:
    	typedef arrayList::value_type value_type;

        arrayList aCenterData;
        arrayList aFaceData[4];
        utPointer utp_data;
        CellData2D();
        CellData2D(int,int=0,int=0,int=0,int=0);
        CellData2D(CellData2D &a);
        CellData2D& operator=(const CellData2D &a);
        void reconstruct(int,int=0,int=0,int=0,int=0);
        bool isEmpty() const;
        ~CellData2D();
        void show_info() const;
    };
    
    typedef CellData2D* pCellData2D;
    
    class CellData3D : public ObjectBase{
    public:
    	typedef arrayList::value_type value_type;

        arrayList aCenterData;
        arrayList aFaceData[6];
        utPointer utp_data;
        CellData3D();
        CellData3D(int,int=0,int=0,int=0,int=0,int=0,int=0);
        CellData3D(CellData3D &a);
        CellData3D& operator=(const CellData3D &a);
        ~CellData3D();
        void show_info() const;
    };
    
    typedef CellData3D* pCellData3D;
    
}

#endif /* CELLDATA_H_ */
