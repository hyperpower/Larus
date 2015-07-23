/*
 * IO_vtk.cpp
 *
 *  Created on: Dec 25, 2014
 *      Author: zhou
 */

#include "IO_vtk.h"
#include "IO.h"
using namespace std;
namespace Larus {

void vtk_unstructured_grid_head(FILE*& data, std::string intru) {
	fprintf(data, "# vtk DataFile Version 3.0\n");
	fprintf(data, "%s \n", intru.c_str());
	fprintf(data, "ASCII\n");
	fprintf(data, "DATASET UNSTRUCTURED_GRID\n");
}

void drawtofile_vtk(std::string filename, Triangle3D& t) {
	FILE* data = open_file(filename, 1);
	vtk_unstructured_grid_head(data, "Triangle3D");
	fprintf(data, "POINTS %d float\n", 3);

	fprintf(data, "%f %f %f\n", t[0].x, t[0].y, t[0].z);
	fprintf(data, "%f %f %f\n", t[1].x, t[1].y, t[1].z);
	fprintf(data, "%f %f %f\n", t[2].x, t[2].y, t[2].z);

	fprintf(data, "\n");
	fprintf(data, "CELLS %d %d \n", 1, 4);
	fprintf(data, "%d %d %d %d\n", 3, 0, 1, 2);

	fprintf(data, "CELL_TYPES %d\n", 1);
	fprintf(data, "%d \n", 5);
	fclose(data);
}

void drawtofile_vtk(std::string filename, Segment3D& t) {
	FILE* data = open_file(filename, 1);
	vtk_unstructured_grid_head(data, "Segment3D");
	fprintf(data, "POINTS %d float\n", 2);

	fprintf(data, "%f %f %f\n", t[0].x, t[0].y, t[0].z);
	fprintf(data, "%f %f %f\n", t[1].x, t[1].y, t[1].z);

	fprintf(data, "\n");
	fprintf(data, "CELLS %d %d \n", 1, 3);
	fprintf(data, "%d %d %d \n", 2, 0, 1);
	fprintf(data, "\n");
	fprintf(data, "CELL_TYPES %d\n", 1);
	fprintf(data, "%d \n", 3);
	fclose(data);
}

void drawtofile_vtk(std::string filename, Point3D& t) {
	FILE* data = open_file(filename, 1);
	vtk_unstructured_grid_head(data, "Segment3D");
	fprintf(data, "POINTS %d float\n", 1);

	fprintf(data, "%f %f %f\n", t.x, t.y, t.z);

	fprintf(data, "\n");
	fprintf(data, "CELLS %d %d \n", 1, 2);
	fprintf(data, "%d %d \n", 1, 0);
	fprintf(data, "\n");
	fprintf(data, "CELL_TYPES %d\n", 1);
	fprintf(data, "%d \n", 1);
	fclose(data);
}

void drawtofile_vtk(std::string filename, ListT<Point3D>& list) {
	FILE* data = open_file(filename, 1);
	vtk_unstructured_grid_head(data, "List Point3D");
	fprintf(data, "POINTS %d float\n", list.size());
	for (ListT<Point3D>::iterator iter = list.begin(); iter != list.end();
			iter++) {
		fprintf(data, "%f %f %f\n", iter->x, iter->y, iter->z);
	}
	fprintf(data, "\n");
	fprintf(data, "CELLS %d %d \n", 1, list.size() + 1);
	fprintf(data, "%d ", list.size());
	for (int i = 0; i < list.size(); i++) {
		fprintf(data, "%d ", i);
	}
	fprintf(data, "\n\n");
	fprintf(data, "CELL_TYPES %d\n", 1);
	fprintf(data, "%d \n", 7); //VTK_POLYGON (=7)
	fclose(data);
}

//===============================================
void vtu_unstructured_grid_file_head(ofstream& fs) {
	fs << "<?xml version=\"1.0\"?>";
	fs
			<< "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">";
	fs << "<UnstructuredGrid>";
}

void vtu_unstructured_grid_file_end(ofstream& fs) {
	fs << "</UnstructuredGrid>";
	fs << "</VTKFile>";
}

void drawtofile_vtu(std::ofstream& fs, ListT<Point3D>& list) {
	for (ListT<Point3D>::iterator iter = list.begin(); iter != list.end();
			iter++) {
		fs << iter->x << " " << iter->y << " " << iter->z << " ";
	}
}
void drawtofile_vtu(std::string filename, Segment3D& t) {
	ofstream fs;
	open_file(fs, filename, 1);
	vtu_unstructured_grid_file_head(fs);
	fs << "<Piece NumberOfPoints=\"" << 2 << "\" NumberOfCells=\"" << 1
			<< "\">";
	fs << "<Points>";
	fs
			<< "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";
	fs << t[0].x << " " << t[0].y << " " << t[0].z << " ";
	fs << t[1].x << " " << t[1].y << " " << t[1].z << " ";
	fs << "</DataArray>";
	fs << "</Points>";
	fs << "<Cells>";
	fs << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
	fs << 0 << " " << 1 << " ";
	fs << "</DataArray>";
	fs << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
	fs << 2;
	fs << "</DataArray>";
	fs << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
	fs << 3;
	fs << "</DataArray>";
	fs << "</Cells>";
	fs << "</Piece>";
	vtu_unstructured_grid_file_end(fs);
	fs.close();
}

}

