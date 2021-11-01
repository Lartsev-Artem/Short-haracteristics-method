#ifndef SHORT_CHARACTERISTICS_BUILD_GRAPH_H
#define SHORT_CHARACTERISTICS_BUILD_GRAPH_H

#include <fstream>
#include <iostream>

#include <set>
#include <string>
#include <vector>

#include <vtk-9.0\vtkCellArray.h>
#include <vtk-9.0\vtkCellData.h>

#include <vtk-9.0\vtkGenericDataObjectReader.h>
#include <vtk-9.0\vtkGenericDataObjectWriter.h>

#include <vtk-9.0\vtkIdList.h>
#include <vtk-9.0/vtkSmartPointer.h>
#include <vtk-9.0\vtkUnstructuredGrid.h>

#include <eigen3/Eigen/Dense>

typedef double Type;
typedef int IntId;
typedef Eigen::Vector3d Vector3;
typedef Eigen::Matrix3d Matrix3;

#include"short_characteristics_calculations.h"  // NormalsToCell

int MainBuildGraphs(int argc, char* argv);

bool CheckCell(const IntId id_cell, const Vector3& direction, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id);

IntId GetNextCell(const Vector3& direction, const std::set<IntId>& boundary_id, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id);

int InitBoundarySet(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::set<IntId>& boundary_cells);
int InitFacesState(const std::vector<IntId>& all_pairs_id, std::vector<IntId>& faces_state);
int ReadGridVtk(const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid);

int ReBuildSetBondary(const IntId id_cell, const Vector3& direction, std::set<IntId>& boundary_id, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id, std::set<IntId>& set_graph);

size_t WriteFileBoundary(const std::string name_file_out, const std::string name_file_graph, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid);


template<typename Str>
size_t ReadStartSettings(Str name_file_settings, Str& name_file_grid, Str& name_file_sphere_direction, Str& name_file_graph) {

	std::ifstream ifile;
	ifile.open(name_file_settings);
	if (!ifile.is_open()) {
		std::cerr << " Error : file settings is not open !\n";
		return 1;
	}

	getline(ifile, name_file_grid);
	getline(ifile, name_file_sphere_direction);
	getline(ifile, name_file_graph);

	ifile.close();
	return 0;
}
#endif