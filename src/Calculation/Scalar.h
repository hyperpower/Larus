/*
 * Scalar.h
 *
 *  Created on: Feb 7, 2015
 *      Author: zhou
 */

#ifndef _SCALAR_H_
#define _SCALAR_H_

#include "../TypeDef.h"
#include "CalDef.h"
#include "../Grid/SPTree.h"
#include "../Grid/SPTreeNode.h"
#include "../Grid/Forest.h"
#include "../Algebra/Arithmetic.h"
#include "../Algebra/Expression.h"

#include <bitset>

namespace Larus
{
//work with tree=================================
//===============================================
//new value======================================
void new_array_on_center_leaf(pOCTree, LarusDef::size_type); //3D
void new_array_on_center_leaf(pQuadTree, LarusDef::size_type); //2D
//shift =========================================
void resize_array_on_center_leaf(pOCTree, LarusDef::size_type); //3D
void resize_array_on_center_leaf(pQuadTree, LarusDef::size_type); //2D
void resize_array_on_center_leaf(Forest2D&, LarusDef::size_type); //2D
void resize_array_on_center_leaf(Forest3D&, LarusDef::size_type); //2D
//-----------------------
void plus_scalar_on_leaf( // 2D tree
		pQuadTree ptree,  //pQuadTree
		arrayList_st& arridx,  //data index
		arrayList& arrval   //data plus
		);
void plus_scalar_on_leaf( // 2D tree
		pOCTree ptree,  //pQuadTree
		arrayList_st& arridx,  //data index
		arrayList& arrval   //data plus
		);
void plus_scalar_on_leaf( // 2D Forest
		Forest2D& forest,  //pQuadTree
		arrayList_st& arridx,  //data index
		arrayList& arrval   //data plus
		);
void plus_scalar_on_leaf(   // 3D Forest
		Forest3D& forest,   //pQuadTree
		arrayList_st& arridx,  //data index
		arrayList& arrval   //data plus
		);
typedef Float (*scalar_pfun)(Float x, Float y, Float z);
void plus_scalar_on_leaf_by_function( // 2D tree
		pQuadTree ptree,  //pQuadTree
		LarusDef::size_type idx,  //data index
		scalar_pfun pf      //data plus
		);
void plus_scalar_on_leaf_by_function( // 2D Forest
		Forest2D& forest,  //pQuadTree
		LarusDef::size_type idx,  //data index
		scalar_pfun pf      //data plus
		);
void set_scalar_on_leaf_by_function( // 2D tree
		pQuadTree ptree,  //pQuadTree
		LarusDef::size_type idx,  //data index
		scalar_pfun pf      //data plus
		);
void set_scalar_on_leaf_by_function( // 2D Forest
		Forest2D& forest,  //pQuadTree
		LarusDef::size_type idx,  //data index
		scalar_pfun pf      //data plus
		);

//set value======================================
//===2D

//===3D========================================
void set_const_on_center_leaf(pOCTree oct, Float f, LarusDef::size_type idx);
void set_const_on_center_leaf(pOCTree oct, arrayList& arr,
		LarusDef::size_type idx_bt);

void set_index_on_center_leaf(Forest2D& forest, LarusDef::size_type = Idx_IDX);
void set_index_on_center_leaf(Forest3D& forest, LarusDef::size_type = Idx_IDX);

void set_value_on_center_leaf(pOCTree, OcTree::pFun_SPTree, utPointer);

typedef Float (*pFun_value)(utPointer, Float, Float, Float);

void set_value_function_on_leaf( //
		Forest2D& forest,    //
		LarusDef::size_type idx, //
		pFun_value pfun,    //
		utPointer utp);

//get value======================================
void getListpNode_leaf_center_data_in_range( // 3D OCTree
		ListT<pOCNode>& listnode,            //as output
		pOCTree pt,                          // ptree
		LarusDef::size_type idx,             //data index
		Float min,                           //min value
		Float max,                           //max value
		TYPE_Range range);                   //range
void getListpNode_leaf_center_data_in_range( // 3D Forest
		ListT<pOCNode>& listnode, //as output
		Forest3D& forest, // Forest
		LarusDef::size_type idx, //data index
		Float min, //min value
		Float max, //max value
		TYPE_Range range);
void getListpNode_leaf_center_data_in_range( // 2D Forest
		ListT<pQTNode>& listnode,           //as output
		Forest2D& forest,                  // Forest
		LarusDef::size_type idx, //data index
		Float min, //min value
		Float max, //max value
		TYPE_Range range);


void getListpNode_leaf_on_line( // 2D Forest
		ListT<pQTNode>& listnode,           //as output
		Forest2D& forest,                  // Forest
		Float loc,
		CSAxis);

void get_average_value_on_level( // 2D QuadTree
		pQuadTree tree,          //Forest
		int level,               //level
		arrayList_st& arridx,    //data index
		ListT<Pair<Point2D, arrayList> >& arrres   //
		);
void get_average_value_on_level( // 2D Forest
		Forest2D& forest,      //Forest
		int level,             //level
		arrayList_st& arridx,  //data index
		ListT<Pair<Point2D, arrayList> >& arrres   //clear the list
		);
void get_max_value(   //
		Forest2D& forest,      //level
		arrayList_st& arridx,  //data index
		arrayList& arrres     //data res
		);
Float get_max_value(   //
		Forest2D& forest,      //level
		LarusDef::size_type idx  //data index
		);
void get_min_value(   //
		Forest2D& forest,      //level
		arrayList_st& arridx,  //data index
		arrayList& arrres     //data res
		);
Float get_min_value(   //
		Forest2D& forest,      //level
		LarusDef::size_type idx  //data index
		);
//interpolate=====================================
void _interpolate_on_axis( // 2D QuadTree Node
		pQTNode pn, //node
		CSAxis axis, //axix
		Float dis, //distance to center of pn
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		);
void _interpolate_1order_on_axis( // 2D QuadTree Node
		pQTNode pn, //node
		CSAxis axis, //axix
		Float dis, //distance to center of pn
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		);
void interpolate_on_face( // 2D QuadTree Node
		pQTNode pn,                         //node
		SPDirection face,                //face
		arrayList_st& arridx, //data index
		arrayList& arrres                   //data res
		);
void interpolate_1order_on_face( // 2D QuadTree Node
		pQTNode pn,                      //node
		SPDirection face,                //face
		arrayList_st& arridx,            //data index
		arrayList& arrres                //data res
		);
void interpolate_1order_on_vertex( // 2D QuadTree Node
		pQTNode pn,                      //node
		SPDirection vertex,              //vertex
		arrayList_st& arridx,            //data index
		arrayList& arrres                //data res
		);
void interpolate_1order_on_vertex_averange_neighbor( // 2D QuadTree Node
		pQTNode pn,                      //node
		SPDirection vertex,              //vertex
		arrayList_st& arridx,            //data index
		arrayList& arrres                //data res
		);
Float interpolate_1order_on_face( // 2D QuadTree Node
		pQTNode pn,                      //node
		SPDirection face,                //face
		LarusDef::size_type idx          //data index
		);
void interpolate_1order_on_face( // 2D QuadTree face
		QTNodeFace& face,                //node
		arrayList_st& arridx,            //data index
		arrayList& arrres                //data res
		);
Float interpolate_1order_on_face( // 2D QuadTree Node
		QTNodeFace& face,                      //face
		LarusDef::size_type           //data index
		);
void interpolate_1order_weight_on_face( // 2D QuadTree Node
		pQTNode pn,                      //node
		SPDirection face,                //face
		arrayList_st& arridx,            //data index
		arrayList& arrres                //data res
		);
void _interpolate_on_vetex( // 2D QuadTree Node
		pQTNode pn, //node
		SPDirection vetex,    //vertex
		arrayList_st& arridx, //data index
		arrayList& arrres   //data res
		);
int _interpolate_node( // 2D QuadTree Node
		pQTNode pn, //node
		const Point2D& point,  //point
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		);
int _interpolate_node_LS( // 2D QuadTree Node Least Square
		pQTNode pn, //node
		const Point2D& point,  //point
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		);
int _interpolate_node_1order( // 2D QuadTree Node
		pQTNode pn, //node
		const Point2D& point,  //point
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		);
int interpolate( // 2D QuadTree
		pQuadTree pqt,  //pQuadTree
		const Point2D& point,  //point
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		);
int interpolate_1order(pQuadTree pqt, const Point2D& point,
		arrayList_st& arridx, arrayList& arrres); //QuadTree
int interpolate(         // 2D Forest
		Forest2D& forest,       //Forest2D
		const Point2D& point,   //point
		arrayList_st& arridx,   //data index
		arrayList& arrres       //data res
		);
int interpolate_1order(
		Forest2D& forest,      //forest2D
		const Point2D& point,  //point
		arrayList_st& arridx,  //idx
		arrayList& arrres);    // 2D Forest
int interpolate_LS(
		Forest2D& forest,      //forest2D
		const Point2D& point,  //point
		arrayList_st& arridx,  //idx
		arrayList& arrres);    // 2D Forest
void interpolate(                  // 2D Forest
		Forest2D& forest,          //Forest
		const ListT<Point2D>& lp,  //point
		arrayList_st& arridx,      //data index
		ListT<Pair<int, arrayList> >& arrres //data res if int=-1 result is wrong
		);
void interpolate(                  //2D Forest
		Forest2D& forest,             //Forest
		const MatrixT<Point2D>& lp,  //point
		arrayList_st& arridx,         //data index
		MatrixT<Pair<int, arrayList> >& arrres //data res if int=-1 result is wrong
		);

//----------
void _interpolate_gradient_on_axis( // 2D QuadTree Node
		pQTNode pn, //node
		CSAxis axis, //axix
		Float dis, //distance to center of pn
		arrayList_st& arridx,  //data index
		arrayList& arrres   //data res
		);
void interpolate_gradient_on_face(    // 2D QuadTree Node
		pQTNode pn,                   //node
		SPDirection face,             //face
		arrayList_st& arridx,         //data index
		arrayList& arrres             //data res
		);
void interpolate_gradient_at_center(  // 2D QuadTree Node
		pQTNode pn,                   //node
		CSAxis axis,                  //axes
		arrayList_st& arridx,         //data index
		arrayList& arrres             //data res
		);
//--------------------------
void interpolate_expression_on_axis( // 2D QuadTree Node
		pQTNode pn,       //node
		CSAxis axis,      //axix
		Float dis,        //distance to center of pn
		Expression& exp   //Expression
		);
void interpolate_expression_on_face( // 2D QuadTree Node
		pQTNode pn,                      //node
		SPDirection face,                //face
		Expression& exp                  //Expression
		);
void interpolate_expression_gradient_on_axis( // 2D QuadTree Node
		pQTNode pn, //node
		CSAxis axis, //axix
		Float dis, //distance to center of pn
		Expression& exp   //Expression
		);
void expression_interpolate_gradient_on_face( // 2D QuadTree Node
		pQTNode pn,                         //node
		SPDirection face,                //face
		Expression& exp   //Expression
		);

//draw===========================================
void draw_vtu_point_leaf_scalars( // 3D OCTree
		std::string filename,   //filename
		pOCTree oct,            //tree
		std::string dataname,   //data name
		LarusDef::size_type idx //data index
		);
void draw_vtu_point_leaf_scalars( // 2D QuadTree
		std::string filename,   //filename
		pQuadTree pqt,          //tree
		std::string dataname,   //data name
		LarusDef::size_type idx //data index
		);
void draw_vtu_point_leaf_scalars( // 2D Forest
		std::string filename,   //filename
		Forest2D& forest,       //forest
		std::string dataname,   //data name
		LarusDef::size_type idx //data index
		);
void draw_vtu_point_leaf_scalars( // 3D OCTree
		std::string filename,  //filename
		pOCTree oct,           //tree
		arrayListT<Pair<std::string, LarusDef::size_type> > arrd  // data pair
		);
void draw_vtu_point_leaf_scalars( // 2D QuadTree
		std::string filename,  //filename
		pQuadTree qt,          //tree
		arrayListT<Pair<std::string, LarusDef::size_type> > arrd  // data pair
		);
void draw_vtu_point_leaf_scalars( // 2D Forest
		std::string filename,  //filename
		Forest2D& forest,      //forest
		arrayListT<Pair<std::string, LarusDef::size_type> > arrd  // data pair
		);
void draw_vtu_as_vector( // 2D QuadTree
		std::string filename,  //filename
		pQuadTree qt,          //tree
		LarusDef::size_type idx_x, //idx x
		LarusDef::size_type idx_y  //icx y
		);
void draw_gnuplot_as_vector( // 2D QuadTree
		std::string filename,  //filename
		pQuadTree qt,          //tree
		int mode,              //mode
		LarusDef::size_type idx_x, //idx x
		LarusDef::size_type idx_y  //icx y
		);
void draw_gnuplot_as_vector(    // 2D Forest
		std::string filename,      //filename
		Forest2D& forest,          //tree
		LarusDef::size_type idx_x, //idx x
		LarusDef::size_type idx_y  //icx y
		);
void draw_gnuplot_as_contour( // 2D QuadTree
		std::string filename,  //filename
		pQuadTree qt,          //tree
		int mode,              //mode
		LarusDef::size_type idx //idx x
		);
void draw_gnuplot_as_contour( // 2D Forest
		std::string filename,   //filename
		Forest2D& forest,       //tree
		int mode,               //mode
		LarusDef::size_type idx //idx x
		);
void draw_gnuplot_as_vector(     // 2D QuadTree
		std::string filename,       //filename
		ListT<pQTNode>& listpnode,  //tree
		int mode,                   //mode
		LarusDef::size_type idx_x,  //idx x
		LarusDef::size_type idx_y   //icx y
		);
void draw_gnuplot_as_contour_line(  // 2D QuadTree
		std::string filename,     //filename
		Forest2D& forest,       //tree
		int mode,               //mode
		LarusDef::size_type idx, //idx x
		arrayList alevel         //array level
		);
void draw_gnuplot_as_contour_line(  // 2D QuadTree
		std::string filename,     //filename
		pQuadTree qt,            //tree
		int mode,               //mode
		LarusDef::size_type idx, //idx x
		arrayList alevel         //array level
		);
void draw_gnuplot_as_contour_line(      // 2D QuadTree
		std::string filename,           //filename
		ListT<pQTNode> lnode,            //tree
		int mode,                       //mode
		LarusDef::size_type idx,        //idx x
		arrayList alevel                //array level
		);
void show_gnuplot_as_contour( // 2D QuadTree
		Forest2D& forest,       //tree
		LarusDef::size_type idx //idx x
		);
void show_gnuplot_as_contour_line( // 2D QuadTree
		Forest2D& forest,       //tree
		LarusDef::size_type idx, //idx x
		LarusDef::size_type nl  //number of level
		);
//check==========================================

}

#endif /* CALULATION_SCALAR_H_ */
