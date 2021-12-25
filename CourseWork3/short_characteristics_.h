#pragma once
using namespace std;

#include <algorithm>
#include <ctime>
#include <execution>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#define PI 3.14159265358979323846
typedef double Type;


#include <vtk-9.0\vtkCellArray.h>
#include <vtk-9.0\vtkCellData.h>
#include <vtk-9.0\vtkDataSet.h>
#include <vtk-9.0\vtkDataSetAttributes.h>
#include <vtk-9.0\vtkDataObject.h>
#include <vtk-9.0\vtkDoubleArray.h>
#include <vtk-9.0\vtkGenericDataObjectReader.h>
#include <vtk-9.0\vtkGenericDataObjectWriter.h>
#include <vtk-9.0\vtkIdList.h>
#include <vtk-9.0\vtkLine.h>
#include <vtk-9.0\vtkLineSource.h>
#include <vtk-9.0\vtkMath.h>
#include <vtk-9.0\vtkNamedColors.h>
#include <vtk-9.0\vtkPointData.h>
#include <vtk-9.0\vtkPoints.h>
#include <vtk-9.0\vtkQuad.h>
#include <vtk-9.0/vtkSmartPointer.h>
#include <vtk-9.0\vtkTetra.h>
#include <vtk-9.0\vtkTriangle.h>
#include <vtk-9.0\vtkUnsignedCharArray.h>
#include <vtk-9.0\vtkUnstructuredGrid.h>

#include <eigen3/Eigen/Dense>






Type CheckLength(const std::vector<std::vector<Type>>& point);
Type Distance(const Type* point1, const Type* point2);
Type Distance(const Type point1_x, const Type point1_y, const Type point1_z,
	const Type point2_x, const Type point2_y, const Type point2_z);

size_t Make2dPoint(const Type* start, Type**& local_basis, const Type* point, Type* new_point);

Type MakeLength(Type* point1, Type* point2);
Type Norm(const Type* vec);
size_t Normalize(Type* vec);
Type Rosh(const Type* S, Type*& a, Type t);
size_t SetBasis(const Type* start_point, const Type* normal, Type* vec_1, Type* vec_2);

size_t SetDirect(const Type* start, const Type* end, Type* direct);
size_t ReadFileVtk(const size_t class_file_vtk, const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, const bool is_print = false);
size_t ReadSphereDirectionVtk(const size_t class_file_vtk, const std::string name_file_sphere_direction, vector<Type>& directions_all);
size_t TransformFileDecartToSphere(const std::string name_file_sphere_direction, const std::string name_new_file_sphere_direction);
size_t TransformNetgenToVtk(const std::string name_file_netgen, const std::string name_new_file_vtk);
size_t SortCellsGrid(Type* main_direction, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vector<int>& sorted_id_cell, vector<Eigen::Vector3d>& centers_tetra);
#include <omp.h>
//#include<tbb/parallel_sort.h>
//#include<tbb/concurrent_vector.h>
size_t PointIntersect(Eigen::Vector3d start_point, Eigen::Vector3d direction, Eigen::Vector3d& point, Eigen::Vector3d& result);