#pragma once
#ifndef SHORT_CHARACTERISTICS_MAIN_H
#define SHORT_CHARACTERISTICS_MAIN_H

#include"short_characteristics_.h"
#include"short_characteristics_calculations.h"
#include "short_characteristics_build_graph.h"
#include<map>

const double eps = 1e-10;

typedef double Type;
typedef Eigen::Vector3d Vector3;
typedef Eigen::Matrix3d Matrix3;

template<typename Type>
size_t ReadStartSettings(std::string name_file_settings, Type& class_file_vtk, std::string& name_file_vtk,
	std::string& name_file_sphere_direction, std::string& out_file_grid_vtk) {

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

	ifile.close();
	return 0;
}

#endif