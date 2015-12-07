/************************
 //  \file   GeoIO_gnuplot.h
 //  \brief
 // 
 //  \author czhou
 //  \date   23 juil. 2014 
 ***********************/
#ifndef IO_GNUPLOT_H_
#define IO_GNUPLOT_H_

#include "../TypeDef.h"
#include <string>
#include "../Geometry/Segment.h"
#include "../Geometry/Polygon.h"
#include "../Geometry/Arrow.h"
#include "../Utility/List.h"
#include "../Utility/ArrayList.h"
#include "../Algebra/MatrixSparCoord.h"
#include "../Grid/Cell.h"
#include "../Grid/Forest.h"
#include "../Calculation/Exp.h"
#include "../Calculation/Poisson.h"
#include "../Calculation/Advection.h"
#include "../Calculation/Scalar.h"

#include "Gnuplot.h"

namespace Larus {

using namespace std;
//==========tmp function
void generate_ListPoint2D(const Point2D& b, const Point2D& e, int num,
		ListT<Point2D>& listp);
void generate_ListPoint2D_circle(const Point2D& c, Float r, int num,
		ListT<Point2D>& listp);

void drawtofile_gnuplot(string filename, const Segment2D& s, int mode);
void drawtofile_gnuplot(string filename, const Point2D& p, int mode);
void drawtofile_gnuplot(string filename, const Arrow2D& p, int mode);
void drawtofile_gnuplot(string filename, const Polygon& p, int mode);

void drawtofile_gnuplot(string filename, const ListT<Point2D>& listp, int mode);
void drawtofile_gnuplot(string filename, const ListT<Segment2D>& listp,
		int mode);
void drawtofile_gnuplot(string filename, const ListT<Arrow2D>&, int mode);
void drawtofile_gnuplot(string filename, const ListT<Pair<Point2D, arrayList> >&, int mode);
void drawtofile_gnuplot(string filename, const arrayListT<Segment2D>&,
		int mode);
void drawtofile_gnuplot(string filename, const arrayListT<Point2D>&, int mode);
void drawtofile_gnuplot(string filename, const Forest2D&, int mode);
//void drawtofile_gnuplot(string filename, ListT<pQTNode> list, int mode);

void drawtofile_gnuplot(string filename, const arrayList& arr, int mode);
void drawtofile_gnuplot(string filename, const arrayList& arr1,
		const arrayList& arr2, int mode);

//void drawtofile_gnuplot(FILE *&data, Point2D p);  //need change

int readfile_gnuplot_countline(string filename);
// ArrayT<Segment2D> readfile_gnuplot(string filename);
// ArrayT<Segment2D> readfile_gnuplot2(string filename) ;

#ifdef MATRIXSPARCOORD_H_
template<class VALUE>
void gnuplot_show(const MatrixSCO<VALUE>& m) {
	Gnuplot gp("boxes");
	arrayList axc(m.NumNonzeros());
	arrayList ayc(m.NumNonzeros());
	arrayList axm(m.NumNonzeros());
	arrayList axp(m.NumNonzeros());
	arrayList aym(m.NumNonzeros());
	arrayList ayp(m.NumNonzeros());
	arrayList val(m.NumNonzeros());
	for (int i = 0; i < m.NumNonzeros(); i++) {
		axc[i] = m.row_ind(i) + 0.5;
		ayc[i] = m.col_ind(i) + 0.5;
		axm[i] = m.row_ind(i);
		axp[i] = m.row_ind(i) + 1;
		aym[i] = m.col_ind(i);
		ayp[i] = m.col_ind(i) + 1;
		val[i] = m.val(i);
	}
	VALUE max = m.max();
	VALUE min = m.min();
	VALUE absmax = ABS(max) > ABS(min) ? ABS(max) : ABS(min);
	string cmdstr = "with boxxy title \"\" fs solid palette";
	gp.set_palette_blue_red();
	gp.set_xrange(0, m.iLen());
	gp.set_yrange_reverse(0, m.jLen());
	gp.set_cbrange(-absmax, absmax);
	std::ostringstream ss;
	ss << "i= " << m.iLen() << " j= " << m.jLen() << " nz= " << m.NumNonzeros();
	gp.set_xlabel(ss.str());
	gp.plot_7(axc, ayc, axm, axp, aym, ayp, val, cmdstr);

}
#endif
template<class VALUE>
void gnuplot_show(const arrayListV<VALUE>& arr) {
	Gnuplot gp("lines");
	gp.plot_x(arr, "");
}
template<class VALUE>
void gnuplot_show(const ListT<VALUE>& list) {
	if (list.size() == 0) {
		return;
	}
	Gnuplot gp("lines");
	arrayListV<VALUE> arr(list.size());
	int i = 0;
	for (typename ListT<VALUE>::const_iterator it = list.begin();
			it != list.end(); ++it) {
		arr[i] = (*it);
		i++;
	}
	gp.plot_x(arr, "");
}
template<class VALUE>
void gnuplot_show_ylog(const ListT<VALUE>& list) {
	if (list.size() == 0) {
		return;
	}
	Gnuplot gp("lines");
	gp.set_ylogscale(10);
	arrayListV<VALUE> arr(list.size());
	int i = 0;
	for (typename ListT<VALUE>::const_iterator it = list.begin();
			it != list.end(); ++it) {
		arr[i] = (*it);
		i++;
	}
	gp.plot_x(arr, "");
}
template<class VALUE>
void gnuplot_show(const MatrixSCR<VALUE>& m) {
	Gnuplot gp("boxes");
	arrayList axc(m.NumNonzeros());
	arrayList ayc(m.NumNonzeros());
	arrayList axm(m.NumNonzeros());
	arrayList axp(m.NumNonzeros());
	arrayList aym(m.NumNonzeros());
	arrayList ayp(m.NumNonzeros());
	arrayList val(m.NumNonzeros());
	int k = 0;
	for (int i = 1; i <= m.iLen(); i++) {
		for (int j = k; j < m.row_ptr(i); j++) {
			axc[k] = i - 1 + 0.5;
			ayc[k] = m.col_ind(k) + 0.5;
			axm[k] = i - 1;
			axp[k] = i;
			aym[k] = m.col_ind(k);
			ayp[k] = m.col_ind(k) + 1;
			val[k] = m.val(k);
			k++;
		}
	}
	VALUE max = m.max();
	VALUE min = m.min();
	VALUE absmax = ABS(max) > ABS(min) ? ABS(max) : ABS(min);
	string cmdstr = "with boxxy title \"\" fs solid palette";
	gp.set_palette_blue_red();
	gp.set_xrange(0, m.iLen());
	gp.set_yrange_reverse(0, m.jLen());
	gp.set_cbrange(-absmax, absmax);
	std::ostringstream ss;
	ss << "i= " << m.iLen() << " j= " << m.jLen() << " nz= " << m.NumNonzeros();
	gp.set_xlabel(ss.str());
	gp.plot_7(axc, ayc, axm, axp, aym, ayp, val, cmdstr);
}

//gnuplot show cell -----------------------------
// gnuplot work with tree -----------------------
inline void gnuplot_inline_data(Gnuplot& gp, const Cell2D& cell) {
	std::ostringstream ss;
	ss << cell.get(CSAxis_X, eCPL_M) << " " << cell.get(CSAxis_Y, eCPL_M)
			<< "\n";
	ss << cell.get(CSAxis_X, eCPL_P) << " " << cell.get(CSAxis_Y, eCPL_M)
			<< "\n";
	ss << cell.get(CSAxis_X, eCPL_P) << " " << cell.get(CSAxis_Y, eCPL_P)
			<< "\n";
	ss << cell.get(CSAxis_X, eCPL_M) << " " << cell.get(CSAxis_Y, eCPL_P)
			<< "\n";
	ss << cell.get(CSAxis_X, eCPL_M) << " " << cell.get(CSAxis_Y, eCPL_M)
			<< "\n";
	gp.cmd(ss.str());
}
void gnuplot_file_data(FILE* data, const Cell2D& cell);

inline void gnuplot_inline_data(Gnuplot& gp, const QTNodeFace& face) {
	std::ostringstream ss;
	switch (face.direction) {
	case SPD_IM: {
		ss << face.pnode->cell->get(CSAxis_X, eCPL_M) << " "
				<< face.pnode->cell->get(CSAxis_Y, eCPL_M) << "\n";
		ss << face.pnode->cell->get(CSAxis_X, eCPL_M) << " "
				<< face.pnode->cell->get(CSAxis_Y, eCPL_P) << "\n";
		break;
	}
	case SPD_IP: {
		ss << face.pnode->cell->get(CSAxis_X, eCPL_P) << " "
				<< face.pnode->cell->get(CSAxis_Y, eCPL_M) << "\n";
		ss << face.pnode->cell->get(CSAxis_X, eCPL_P) << " "
				<< face.pnode->cell->get(CSAxis_Y, eCPL_P) << "\n";
		break;
	}
	case SPD_JP: {
		ss << face.pnode->cell->get(CSAxis_X, eCPL_M) << " "
				<< face.pnode->cell->get(CSAxis_Y, eCPL_P) << "\n";
		ss << face.pnode->cell->get(CSAxis_X, eCPL_P) << " "
				<< face.pnode->cell->get(CSAxis_Y, eCPL_P) << "\n";
		break;
	}
	case SPD_JM: {
		ss << face.pnode->cell->get(CSAxis_X, eCPL_M) << " "
				<< face.pnode->cell->get(CSAxis_Y, eCPL_M) << "\n";
		ss << face.pnode->cell->get(CSAxis_X, eCPL_P) << " "
				<< face.pnode->cell->get(CSAxis_Y, eCPL_M) << "\n";
		break;
	}
	default:{
		ASSERT_MSG(false, " !> Error face direction (gnuplot_inline_data)\n");
	}
	}
	gp.cmd(ss.str());
}
inline void gnuplot_show(const pQTNode pnode) {
	Gnuplot gp("lines");
	std::ostringstream ss;
	ss << "plot \"-\" using 1:2 title \"\" " << "with lines lw 1";
	gp.set_equal_ratio();
	gp.cmd(ss.str());
	ss.str("");
	gnuplot_inline_data(gp, *(pnode->cell));
	gp.cmd("e");
}
inline void gnuplot_show(const ListT<pQTNode>& lpnode) {
	Gnuplot gp("lines");
	std::ostringstream ss;
	ss << "plot \"-\" using 1:2 title \"\" " << "with lines lw 1";
	gp.set_equal_ratio();
	gp.cmd(ss.str());
	ss.str("");
	for (typename ListT<pQTNode>::const_iterator it = lpnode.begin();
			it != lpnode.end(); ++it) {
		gnuplot_inline_data(gp, *((*it)->cell));
		gp.cmd("\n");
	}
	gp.cmd("e");
}
inline void gnuplot_show(const Forest2D& forest) {
	Gnuplot gp("lines");
	std::ostringstream ss;
	ss << "plot \"-\" using 1:2 title \"\" " << "with lines lw 1";
	gp.set_equal_ratio();
	gp.cmd(ss.str());
	ss.str("");
	for (typename Forest2D::const_iterator it = forest.begin();
			it != forest.end(); ++it) {
		gnuplot_inline_data(gp, (*it->cell));
		gp.cmd("\n");
	}
	gp.cmd("e");
}
inline void gnuplot_show_different_level(const Forest2D& forest) {
	Gnuplot gp("lines");
	std::ostringstream ss;
	ss << "plot \"-\" using 1:2 title \"\" " << "with lines lw 1";
	gp.set_equal_ratio();
	gp.set_xrange(6, 6.5);
	gp.cmd(ss.str());
	ss.str("");
	for (typename Forest2D::const_iterator_face it = forest.begin_face();
			it != forest.end_face(); ++it) {
		if (it->face_type != SPFT_Equal) {
			gnuplot_inline_data(gp, (*it));
			gp.cmd("\n");
		}
	}
	gp.cmd("e");
}

void gnuplot_file_data(FILE* data, const QTNodeFace& face);

inline void gnuplot_datafile_different_level(string filename, const Forest2D& forest) {
	FILE* data = open_file(filename, 1);
	for (typename Forest2D::const_iterator_face it = forest.begin_face();
			it != forest.end_face(); ++it) {
		if (it->face_type != SPFT_Equal) {
			gnuplot_file_data(data , (*it));
			fprintf(data, "\n");
		}
	}
	fclose(data);
}

//------------------------------------------------
inline void gnuplot_show(const Exp2D& exp) {
	Gnuplot gp("lines");
	for (typename Exp2D::const_iterator it = exp.begin(); it != exp.end();
			++it) {
		Point2D cp = it->pnode->cell->getCenterPoint();
		std::ostringstream ss;
		ss << "\"" << " ( " << it->pnode->data->aCenterData[Idx_IDX] << ", "
				<< parseQTNodeType(it->pnode->getType()) << " )\" at first "
				<< cp.x << ", first " << cp.y << " center";
		gp.set_label(ss.str());

	}
	std::ostringstream ss;
	ss << "plot \"-\" using 1:2 title \"\" " << "with lines lw 1";
	gp.set_equal_ratio();
	gp.cmd(ss.str());
	ss.str("");
	for (typename Exp2D::const_iterator it = exp.begin(); it != exp.end();
			++it) {
		if (it->idx != ExpTerm::IDX_CONST) {
			gnuplot_inline_data(gp, (*it->pnode->cell));
			gp.cmd("\n");
		}
	}
	gp.cmd("e");

}

inline void gnuplot_show_as_contour(const Forest2D& forest, int cv_idx) {
	typedef Forest2D::Data::value_type vt;
	ListT<vt> lxc, lyc, lxm, lxp, lym, lyp, lval;
	for (typename Forest2D::const_iterator iter = forest.begin();
			iter != forest.end(); iter++) {
		const typename Forest2D::Node* pnode = iter.get_pointer();
		lxc.push_back(pnode->cell->getCenterPoint().x);
		lyc.push_back(pnode->cell->getCenterPoint().y);
		lxm.push_back(pnode->cell->getMM().x);
		lxp.push_back(pnode->cell->getPP().x);
		lym.push_back(pnode->cell->getMM().y);
		lyp.push_back(pnode->cell->getPP().y);
		lval.push_back(pnode->data->aCenterData[cv_idx]);
	}
	typename ListT<vt>::const_iterator iter = lval.begin();
	vt max = (*iter);
	vt min = (*iter);
	for (++iter; iter != lval.end(); ++iter) {
		if ((*iter) > max) {
			max = (*iter);
		}
		if ((*iter) < min) {
			min = (*iter);
		}
	}
	Gnuplot gp("boxes");
	string cmdstr = "with boxxy title \"\" fs solid palette";
	gp.set_palette_blue_red();
	if (max == min) {
		gp.set_cbrange(-1, 1);
	} else {
		gp.set_cbrange(min, max);
	}
	gp.set("format cb \"%+-.2e\"");
	std::ostringstream ss;
	gp.set_xlabel(ss.str());
	gp.set_equal_ratio();
	gp.plot_7(lxc, lyc, lxm, lxp, lym, lyp, lval, cmdstr);
}

inline void gnuplot_show_as_surface(Forest2D& forest, int cv_idx, string title =
		"") {
	typedef Forest2D::Data::value_type vt;
	typedef Forest2D::Cell Cell;
	typedef Forest2D::Node Node;
	Gnuplot gp("lines");
	std::ostringstream ss;
	ss << "splot \"-\" using 1:2:3 title \"\" " << "with lines lw 1 lc palette";
	gp.set_palette_blue_red();
	gp.set_title(title);
	// x rot , z rot , 1 scale, 1 zscale
	gp.set("view 40,10, 1.1, 1");
	gp.set("ticslevel 0");
	gp.set("colorbox");
	gp.set_equal_ratio();
	gp.cmd(ss.str());
	ss.str("");
	for (typename Forest2D::iterator it = forest.begin(); it != forest.end();
			++it) {
		Node* pn = it.get_pointer();
		arrayList aval(4), ares(1);
		arrayList_st aidx(1);
		aidx[0] = cv_idx;
		const SPDirection index_dir[4] = { SPD_MM, SPD_PM, SPD_PP, SPD_MP };
		//  0       1       2       3
		arrayListT<Point2D> aloc(4);
		for (int i = 0; i < 4; ++i) {
			//get value on vertex =====================
			interpolate_1order_on_vertex_averange_neighbor(pn, index_dir[i],
					aidx, ares);
			aval[i] = ares[0];
			aloc[i] = pn->getPoint(index_dir[i]); //get point location
		}
		std::ostringstream ss;
		for (int i = 0; i < 4; ++i) {
			ss << aloc[i].x << " " << aloc[i].y << " " << aval[i] << "\n";
		}
		ss << aloc[0].x << " " << aloc[0].y << " " << aval[0] << "\n";
		gp.cmd(ss.str());
		gp.cmd("\n");
	}
	gp.cmd("e");
}

#ifdef _POISSON_H_
inline void gnuplot_show_as_contour(const Poisson_Eq<Dimension_2D>& pe) {
	gnuplot_show_as_contour((*pe.pforest), pe.phi_idx);
}
#endif

#ifdef _ADVECTION_H_
inline void gnuplot_show_as_contour(const Advection_Eq<Dimension_2D>& ae) {
	gnuplot_show_as_contour((*ae.pforest), ae.phi_idx);
}
#endif

//-----------------------------------------------

}//end namespace

#endif /* GEOIO_GNUPLOT_H_ */
