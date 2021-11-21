#pragma once
#ifndef SHORT_CHARACTERISTICS_MAIN_H
#define SHORT_CHARACTERISTICS_MAIN_H

#include"short_characteristics_.h"
#include"short_characteristics_calculations.h"
#include "short_characteristics_build_graph.h"
#include<map>

int ReadGraph(const std::string name_file_graph, std::vector<IntId>& sorted_id_cell);
int Start(const std::string name_file_graph, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid,
	const std::vector<int>& all_pairs_face, std::map<int, Vector3>& nodes_value,
	const vector<Type>& directions, const vector<Type>& squares, const Type square_surface,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate,
	std::vector<Type>& Illum1);


const double eps = 1e-10;

typedef double Type;
typedef Eigen::Vector3d Vector3;
typedef Eigen::Matrix3d Matrix3;

template<typename Type>
size_t ReadStartSettings(std::string name_file_settings, Type& class_file_vtk, std::string& name_file_vtk,
	std::string& name_file_sphere_direction, std::string& out_file_grid_vtk, std::string& name_file_graph, std::string& out_file_E1d) {

	std::ifstream ifile;
	ifile.open(name_file_settings);
	if (!ifile.is_open()) {
		std::cerr << " Error : file settings build graph is not open !\n";
		return 1;
	}

	std::string str; // переменная для перевода строки при чтении из файла

	ifile >> class_file_vtk;
	getline(ifile, str);
	getline(ifile, name_file_vtk);
	getline(ifile, name_file_sphere_direction);
	getline(ifile, out_file_grid_vtk);
	getline(ifile, name_file_graph);
	getline(ifile, out_file_E1d);

	ifile.close();
	return 0;
}

#endif