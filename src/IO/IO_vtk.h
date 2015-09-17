/*
 * IO_vtk.h
 *
 *  Created on: Dec 25, 2014
 *      Author: zhou
 */

#ifndef IO_VTK_H_
#define IO_VTK_H_

#include "../TypeDef.h"
#include <string>
#include <fstream>
#include "../Geometry/Polygon.h"
#include "../Utility/List.h"
#include "../Utility/ArrayList.h"

#include "../Geometry/Triangle.h"




namespace Larus {

void vtk_unstructured_grid_head(FILE*& data, std::string intru);

void vtu_unstructured_grid_file_head(std::ofstream& fs);
void vtu_unstructured_grid_file_end(std::ofstream& fs);

void drawtofile_vtu(std::ofstream& fs, ListT<Point3D>&);

void drawtofile_vtk(std::string filename, Triangle3D&);
void drawtofile_vtk(std::string filename, Segment3D&);
void drawtofile_vtk(std::string filename, Point3D&);
void drawtofile_vtk(std::string filename, ListT<Point3D>&);

void drawtofile_vtu(std::string filename, Segment3D&);

void drawtofile_vtu(std::ofstream& fs, ListT<Point3D>&);
/*
 * This part is based on vtk library > 5.8
 */
int vtk_show(const Segment3D& seg);

}

#endif /* IO_IO_VTK_H_ */
