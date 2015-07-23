/************************
 //  \file   Point.h
 //  \brief
 // 
 //  \author czhou
 //  \date   21 mai 2015 
 ***********************/
#ifndef TS_POINT_H_
#define TS_POINT_H_

#include <array>
#include <assert.h>
#include "ts_define.h"
#include "ts_predicates.h"

namespace LarusTS {

//#define     GTS_POINT_CLASS                 (klass)
//#define     GTS_POINT                       (obj)
//#define     GTS_IS_POINT                    (obj)
//            GtsPointClass;
//            GtsPoint;
//
//GtsPointClass* gts_point_class              (void);
//GtsPoint*   gts_point_new                   (GtsPointClass *klass,
//                                             gdouble x,
//                                             gdouble y,
//                                             gdouble z);
//void        gts_point_set                   (GtsPoint *p,
//                                             gdouble x,
//                                             gdouble y,
//                                             gdouble z);
//#define     gts_point_is_in_rectangle       (p, p1, p2)
//GtsPoint*   gts_segment_triangle_intersection
//                                            (GtsSegment *s,
//                                             GtsTriangle *t,
//                                             gboolean boundary,
//                                             GtsPointClass *klass);
//void        gts_point_transform             (GtsPoint *p,
//                                             GtsMatrix *m);
//gdouble     gts_point_distance              (GtsPoint *p1,
//                                             GtsPoint *p2);
//gdouble     gts_point_distance2             (GtsPoint *p1,
//                                             GtsPoint *p2);
//gdouble     gts_point_orientation_3d        (GtsPoint *p1,
//                                             GtsPoint *p2,
//                                             GtsPoint *p3,
//                                             GtsPoint *p4);
//gint        gts_point_orientation_3d_sos    (GtsPoint *p1,
//                                             GtsPoint *p2,
//                                             GtsPoint *p3,
//                                             GtsPoint *p4);
//enum        GtsIntersect;
//gdouble     gts_point_in_circle             (GtsPoint *p,
//                                             GtsPoint *p1,
//                                             GtsPoint *p2,
//                                             GtsPoint *p3);
//gdouble     gts_point_in_triangle_circle    (GtsPoint *p,
//                                             GtsTriangle *t);
//GtsIntersect gts_point_is_in_triangle       (GtsPoint *p,
//                                             GtsTriangle *t);
//gdouble     gts_point_orientation           (GtsPoint *p1,
//                                             GtsPoint *p2,
//                                             GtsPoint *p3);
//gint        gts_point_orientation_sos       (GtsPoint *p1,
//                                             GtsPoint *p2,
//                                             GtsPoint *p3);
//gdouble     gts_point_segment_distance2     (GtsPoint *p,
//                                             GtsSegment *s);
//gdouble     gts_point_segment_distance      (GtsPoint *p,
//                                             GtsSegment *s);
//void        gts_point_segment_closest       (GtsPoint *p,
//                                             GtsSegment *s,
//                                             GtsPoint *closest);
//gdouble     gts_point_triangle_distance     (GtsPoint *p,
//                                             GtsTriangle *t);
//void        gts_point_triangle_closest      (GtsPoint *p,
//                                             GtsTriangle *t,
//                                             GtsPoint *closest);
//gdouble     gts_point_triangle_distance2    (GtsPoint *p,
//                                             GtsTriangle *t);
//gboolean    gts_point_is_inside_surface     (GtsPoint *p,
//                                            GNode *tree,
//                                            gboolean is_open);

template<class TYPE, std::size_t DIM> class Segment;
template<class TYPE, std::size_t DIM> class Surface;
template<class TYPE, std::size_t DIM> class Edge;
template<class TYPE, std::size_t DIM> class Vertex;
template<class TYPE, std::size_t DIM> class Triangle;

template<class TYPE, std::size_t DIM>
class Point: public std::array<TYPE, DIM> {
public:
	typedef std::array<TYPE, DIM> base_class;
	typedef TYPE value_type;
	typedef value_type* pointer;
	typedef const value_type* const_pointer;
	typedef value_type& reference;
	typedef const value_type& const_reference;
	typedef value_type* iterator;
	typedef const value_type* const_iterator;
	typedef std::size_t size_type;
	typedef std::ptrdiff_t difference_type;
	typedef std::reverse_iterator<iterator> reverse_iterator;
	typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

	const static int Dim = DIM;

	//constructor
	Point() :
		base_class() {
	}

	Point(const TYPE& a, const TYPE& b, const TYPE& c = 0) :
		base_class() {
		this->__elems_[0] = a;
		this->__elems_[1] = b;
		if (DIM == 3) {
			this->__elems_[2] = c;
		}
	}

	void set(const TYPE& a, const TYPE& b, const TYPE& c = 0) {
		this->__elems_[0] = a;
		this->__elems_[1] = b;
		if (DIM == 3) {
			this->__elems_[2] = c;
		}
	}

	const_reference x() const {
		return this->__elems_[0];
	}

	reference x() {
		return this->__elems_[0];
	}

	const_reference y() const {
		return this->__elems_[1];
	}

	reference y() {
		return this->__elems_[1];
	}

	const_reference z() const {
		assert(DIM == 3);
		return this->__elems_[2];
	}

	reference z() {
		assert(DIM == 3);
		return this->__elems_[2];
	}

	void reconstruct(const TYPE& a, const TYPE& b, const TYPE& c = 0) {
		this->__elems_[0] = a;
		this->__elems_[1] = b;
		if (DIM == 3) {
			this->__elems_[2] = c;
		}
	}

	bool operator==(const Point<TYPE, DIM> &a) const {
		if (DIM == 2) {
			return (this->__elems_[0] == a[0] && this->__elems_[1] == a[1]) ?
					true : false;
		} else {
			return (this->__elems_[0] == a[0] && this->__elems_[1] == a[1]
					&& this->__elems_[2] == a[2]) ? true : false;
		}
	}
	bool operator!=(const Point<TYPE, DIM> &a) const {
		if (DIM == 2) {
			return !(
					(this->__elems_[0] == a[0] && this->__elems_[1] == a[1]) ?
							true : false);
		} else {
			return !(
					(this->__elems_[0] == a[0] && this->__elems_[1] == a[1]
							&& this->elems[2] == a[2]) ? true : false);
		}
	}
	void show() const {
		std::cout << std::scientific << "( " << this->__elems_[0] << " , "
				<< this->__elems_[1];
		if (DIM == 3) {
			std::cout << " , " << this->__elems_[2] << " )\n";
		} else {
			std::cout << " )\n";
		}
	}

	template<typename T>
	void transfer(const T&dx, const T&dy, const T&dz) {
		this->__elems_[0] = this->__elems_[0] + TYPE(dx);
		this->__elems_[1] = this->__elems_[1] + TYPE(dy);
		if (DIM == 3) {
			this->__elems_[2] = this->__elems_[2] + TYPE(dz);
		}
	}

	template<typename T>
	void scale(const T&dx, const T&dy, const T&dz) {
		this->__elems_[0] = this->__elems_[0] * TYPE(dx);
		this->__elems_[1] = this->__elems_[1] * TYPE(dy);
		if (DIM == 3) {
			this->__elems_[2] = this->__elems_[2] * TYPE(dz);
		}
	}

	inline size_type size() const {
		return size_type(DIM);
	}

};

//function out of class
template<class POINT>
typename POINT::value_type point_distance(const POINT& p1, const POINT& p2) {
	if (POINT::Dim == 2) {
		return sqrt(
				(p1[0] - p2[0]) * (p1[0] - p2[0])
						+ (p1[1] - p2[1]) * (p1[1] - p2[1]));
	} else {
		return sqrt(
				(p1[0] - p2[0]) * (p1[0] - p2[0])
						+ (p1[1] - p2[1]) * (p1[1] - p2[1])
						+ (p1[2] - p2[2]) * (p1[2] - p2[2]));
	}
}
template<class POINT>
typename POINT::value_type point_distance2(const POINT& p1, const POINT& p2) {
	if (POINT::Dim == 2) {
		return (p1[0] - p2[0]) * (p1[0] - p2[0])
				+ (p1[1] - p2[1]) * (p1[1] - p2[1]);
	} else {
		return (p1[0] - p2[0]) * (p1[0] - p2[0])
				+ (p1[1] - p2[1]) * (p1[1] - p2[1])
				+ (p1[2] - p2[2]) * (p1[2] - p2[2]);
	}
}
template<class POINT>
bool point_is_in_on_rectangle(const POINT& p, const POINT& p1,
		const POINT& p2) {
	bool res = true;
	for (typename POINT::size_type i = 0; i < POINT::Dim; ++i) {
		res = res && p[i] >= p1[i] && p[i] <= p2[i];
	}
}
/**
 * gts_point_orientation:
 * @p1: a #GtsPoint.
 * @p2: a #GtsPoint.
 * @p3: a #GtsPoint.
 *
 * Checks for orientation of the projection of three points on the
 * (x,y) plane. The result is also an approximation of twice the
 * signed area of the triangle defined by the three points. This
 * function uses adaptive floating point arithmetic and is
 * consequently geometrically robust.
 *
 * Returns: a positive value if @p1, @p2 and @p3 appear in
 * counterclockwise order, a negative value if they appear in
 * clockwise order and zero if they are colinear.
 */
template<class POINT>
Float point_orientation(const POINT& p1, const POINT& p2, const POINT& p3) {
	return orient2d((double *) p1.data(), (double *) p2.data(),
			(double *) p3.data());
}

/**
 * point_orientation_3d:
 * @p1: a Point.
 * @p2: a Point.
 * @p3: a Point.
 * @p4: a Point.
 *
 * Checks if @p4 lies above, below or on the plane passing through the
 * points @p1, @p2 and @p3. Below is defined so that @p1, @p2 and @p3
 * appear in counterclockwise order when viewed from above the
 * plane. The returned value is an approximation of six times the
 * signed volume of the tetrahedron defined by the four points. This
 * function uses adaptive floating point arithmetic and is
 * consequently geometrically robust.
 *
 * Returns:
 *        a positive value if @p4 lies below,
 *        a negative value if @p4 lies above the plane,
 *        zero             if the four points are coplanar.
 */
template<class POINT>
Float point_orientation_3d(const POINT& p1, const POINT& p2, const POINT& p3,
		const POINT& p4) {
	assert(POINT::Dim == 3);
	return orient3d((double *) p1.data(), (double *) p2.data(),
			(double *) p3.data(), (double *) p4.data());
}
/**
 * point_is_in_triangle:
 * @p: a #GtsPoint.
 * @t: a #GtsTriangle.
 *
 * Tests if the planar projection (x, y) of @p is inside, outside or
 * on the boundary of the planar projection of @t.  This function is
 * geometrically robust.
 *
 * Returns: %GTS_IN if @p is inside @t, %GTS_ON if @p is on the boundary of
 * @t, %GTS_OUT otherwise.
 */
template<class TYPE, std::size_t DIM>
Intersect point_is_in_triangle(const Point<TYPE, DIM>& p,
		const Triangle<TYPE, DIM>& t) {
	Vertex<TYPE, DIM> v1, v2, v3;
	Float d1, d2, d3;

	triangle_vertices(t, &v1, &v2, &v3);  //

	d1 = point_orientation((v1), (v2), p);
	if (d1 < 0.0)
		return TS_OUT;
	d2 = point_orientation((v2), (v3), p);
	if (d2 < 0.0)
		return TS_OUT;
	d3 = point_orientation((v3), (v1), p);
	if (d3 < 0.0)
		return TS_OUT;
	if (d1 == 0.0 || d2 == 0.0 || d3 == 0.0)
		return TS_ON;
	return TS_IN;
}

/**
 * gts_point_in_triangle_circle:
 * @p: a #GtsPoint.
 * @t: a #GtsTriangle.
 *
 * Tests if the planar projection (x, y) of @p is inside or outside
 * the circumcircle of the planar projection of @t. This function is
 * geometrically robust.
 *
 * Returns: a positive number if @p lies inside,
 * a negative number if @p lies outside and zero if @p lies on
 * the circumcircle of @t.
 */
template<class TYPE, std::size_t DIM>
Float point_in_triangle_circle(     //
		const Point<TYPE, DIM>& p,  //
		const Triangle<TYPE, DIM>& t) {
	Vertex<TYPE, DIM> v1, v2, v3;
	triangle_vertices(t, v1, v2, v3);
	return incircle((double *) v1.data(), (double *) v2.data(),
			(double *) v3.data(), (double *) &p.data());
}

/**
 * point_in_circle:
 * @p: a #GtsPoint.
 * @p1: a #GtsPoint.
 * @p2: a #GtsPoint.
 * @p3: a #GtsPoint.
 *
 * Tests if the planar projection (x, y) of @p is inside or outside the
 * circle defined by the planar projection of @p1, @p2 and @p3.
 *
 * Returns: a positive number if @p lies inside,
 * a negative number if @p lies outside and zero if @p lies on
 * the circle.
 */
template<class TYPE, std::size_t DIM>
Float gts_point_in_circle(const Point<TYPE, DIM>& p,   //
		const Point<TYPE, DIM>& p1,  //
		const Point<TYPE, DIM>& p2,  //
		const Point<TYPE, DIM>& p3) {
	return incircle((double *) p1.data(), (double *) p2.data(),
			(double *) p3.data(), (double *) p.data());
}
template<class TYPE, std::size_t DIM>
Float gts_point_in_sphere(  //
		const Point<TYPE, DIM>& p,  //
		const Point<TYPE, DIM>& p1, //
		const Point<TYPE, DIM>& p2,  //
		const Point<TYPE, DIM>& p3, //
		const Point<TYPE, DIM>& p4) {
	return insphere((double *) &p1.data(), (double *) &p2.data(),
			(double *) &p3.data(), (double *) &p4.data(), (double *) &p.data());
}

/**
 * gts_point_segment_distance2:
 * @p: a #GtsPoint.
 * @s: a #GtsSegment.
 *
 * Returns: the square of the minimun Euclidean distance between @p and @s.
 */
template<class TYPE, std::size_t DIM>
TYPE point_segment_distance2( //
		const Point<TYPE, DIM>& p, //
		const Segment<TYPE, DIM>& s) {
	TYPE t, ns2, x, y, z;
	const Point<TYPE, DIM>* p1, p2;

	p1 = s->v1;
	p2 = s->v2;
	ns2 = point_distance2(*p1, *p2);
	if (ns2 == 0.0)
		return point_distance2(*p, *p1);
	t = ((p2->x() - p1->x()) * (p->x() - p1->x())
			+ (p2->y() - p1->y()) * (p->y() - p1->y())
			+ (p2->z() - p1->z()) * (p->z() - p1->z())) / ns2;
	if (t > 1.0)
		return point_distance2(*p, *p2);
	if (t < 0.0)
		return point_distance2(*p, *p1);
	x = (1. - t) * p1->x() + t * p2->x() - p->x();
	y = (1. - t) * p1->y() + t * p2->y() - p->y();
	z = (1. - t) * p1->z() + t * p2->z() - p->z();
	return x * x + y * y + z * z;
}

/**
 * gts_point_segment_distance:
 * @p: a #GtsPoint.
 * @s: a #GtsSegment.
 *
 * Returns: the minimun Euclidean distance between @p and @s.
 */
template<class TYPE, std::size_t DIM>
Float point_segment_distance( //
		const Point<TYPE, DIM>& p, //
		const Segment<TYPE, DIM>& s) {
	return sqrt(point_segment_distance2(p, s));
}

/**
 * gts_point_segment_closest:
 * @p: a #GtsPoint.
 * @s: a #GtsSegment.
 * @closest: a #GtsPoint.
 *
 * Set the coordinates of @closest to the coordinates of the point belonging
 * to @s closest to @p.
 */
template<class TYPE, std::size_t DIM>
void gts_point_segment_closest(const Point<TYPE, DIM>& p,
		const Segment<TYPE, DIM>& s, Point<TYPE, DIM>& closest) {
	TYPE t, ns2;
	Point<TYPE, DIM>* p1, *p2;

	p1 = s->v1;
	p2 = s->v2;
	ns2 = point_distance2(p1, p2);

	if (ns2 == 0.0) {
		closest = (*p1);
		return;
	}

	t = ((p2->x() - p1->x()) * (p->x() - p1->x())
			+ (p2->y() - p1->y()) * (p->y() - p1->y())
			+ (p2->z() - p1->z()) * (p->z() - p1->z())) / ns2;

	if (t > 1.0)
		closest = (*p2);
	else if (t < 0.0)
		closest = (*p1);
	else
		closest.set( //
				(1. - t) * p1->x() + t * p2->x(), //
				(1. - t) * p1->y() + t * p2->y(), //
				(1. - t) * p1->z() + t * p2->z());
}

/**
 * gts_point_orientation_sos:
 * @p1: a #GtsPoint.
 * @p2: a #GtsPoint.
 * @p3: a #GtsPoint.
 *
 * Checks for orientation of the projection of three points on the
 * (x,y) plane.
 *
 * Simulation of Simplicity (SoS) is used to break ties when the
 * orientation is degenerate (i.e. @p3 lies on the line defined by
 * @p1 and @p2).
 *
 * Returns: a positive value if @p1, @p2 and @p3 appear in
 * counterclockwise order or a negative value if they appear in
 * clockwise order.
 */
template<class TYPE, std::size_t DIM>
int point_orientation_sos( //
		const Point<TYPE, DIM>& p1, //
		const Point<TYPE, DIM>& p2,  //
		const Point<TYPE, DIM>& p3 //
		) {
	Float o;

	o = orient2d((double *) p1.data(), (double *) p2.data(), (double *) p3.data());
	if (o != 0.)
		return SIGN(o);
	else {
		const Point<TYPE, DIM>* p[3];
		int sign;

		p[0] = &p1;
		p[1] = &p2;
		p[2] = &p3;

		sign = sortp(p, 3);

		/* epsilon^1/4 */
		o = ORIENT1D(p[1]->x(), p[2]->x());
		if (o != 0.)
			return -SIGN(o) * sign;

		/* epsilon^1/2 */
		o = ORIENT1D(p[1]->y(), p[2]->y());
		if (o != 0.)
			return SIGN(o) * sign;

		/* epsilon */
		o = ORIENT1D(p[0]->x(), p[2]->x());
		if (o != 0.)
			return SIGN(o) * sign;

		/* epsilon^3/2 */
		return sign;
	}
}

} //end of namespace

#endif /* POINT_H_ */
