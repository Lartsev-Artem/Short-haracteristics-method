#pragma once
#ifndef SHORT_CHARACTERISTICS_CALCULATIONS_H
#define SHORT_CHARACTERISTICS_CALCULATIONS_H

#include"short_characteristics_.h"
#include "short_characteristics_main.h"

#include<map>
typedef double Type;
typedef Eigen::Vector3d Vector3;
typedef Eigen::Matrix3d Matrix3;
int InitGlobalValue(Vector3& start_point_plane_coord, Matrix3& transform_matrix, Matrix3& inverse_transform_matrix,
	Matrix3& straight_face, Matrix3& inclined_face);

Type CalculateIllumeOnInnerFace(const int num_in_face, const int global_num_in_face, const Eigen::Matrix4d& vertex_tetra,
	const Vector3& x, const Vector3& x0, const std::map<int, Vector3>& nodes_value,
	const Matrix3& straight_face, const Matrix3& inclined_face, const Matrix3& transform_matrix, const Vector3& start_point_plane_coord);

size_t SortCellsGrid(Vector3 main_direction, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vector<int>& sorted_id_cell, 
	const vector<Vector3>& centers_tetra);
Type IntegarteDirection(const vector<Type>& directions);
Type BoundaryFunction(const Vector3 x);
int ReadNextCellData(ifstream& ifs, Eigen::Matrix<Type, 4, 3>& normals, Vector3& center);
int ReadCenters(const std::string name_file, vector<Vector3>& centers_tetra);
int ReadCentersFromGeneralFile(const std::string name_file, vector<Vector3>& centers_tetra);

Type GetValueInCenterCell(const int num_cell, vtkCell* cur_cell, const Vector3 center, const Vector3 direction,
	const Eigen::Matrix4d& vertex_tetra, const std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate,
	const Matrix3& straight_face, const Matrix3& inclined_face, const Matrix3& transform_matrix, const Vector3& start_point_plane_coord);

int FromTetraToPlane(const Eigen::Matrix3d& transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& tetra_coord, Eigen::Vector3d& plane_coord);

int GetNodesValues(const int num_cur_cell, vtkCell* cur_cell, const int num_cur_out_face, const int* face_state,
	const Vector3& direction, const Eigen::Matrix4d& vertex_tetra,
	std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate,
	const Matrix3& straight_face, const Matrix3& inclined_face, const Matrix3& transform_matrix, const Matrix3& inverse_transform_matrix,
	const Vector3& start_point_plane_coord);


int SetVertexMatrix(vtkCell* cell, Eigen::Matrix4d& vertex_tetra);
int FindInAndOutFaces(const Vector3  direction, vtkCell* cur_cell_tetra, int* face_state);
int FindNeighborsPairFace(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<int>& all_pairs_face);
int SetVertexMatrix(const size_t number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, Eigen::Matrix4d& vertex_tetra);
int InitNodesValue(const std::vector<int>& all_pairs_face, std::map<int, Vector3>& nodes_value);
int MakeFileCellData(const std::string name_file, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, bool rewrite_if_exist = true);
int GetNumberNeighborFace(const int a, const int b, const int c, vtkCell* neighbor_cell);
int MakeFileDirections(const std::string name_file_sphere_direction, vector<Type>& directions_all);
int SetVertexMatrix(ifstream& ifs, const size_t number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, Eigen::Matrix4d& vertex_tetra);
int FromLocalToGlobalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& local_coord, Eigen::Vector3d& global_coord);
int FromGlobalToLocalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& global_coord, Eigen::Vector3d& local_coord);

int GetInterpolationCoef    (const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value, Eigen::Vector3d& coef);
Vector3 GetInterpolationCoef(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value);
Type GetIllum(const int cur_id, const Vector3 x, const Type s, const Type I_node_prev,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate);
int NormalToFace(vtkCell* cell_face, Vector3& n);
size_t IntersectionWithPlane(vtkCell* face, const Vector3& start_point, const Vector3& direction, Vector3& result);
bool InTriangle(vtkCell* face, const Eigen::Vector3d& XX);
int FindCentersSpheres(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, vector<Vector3>& centers_tetra);
int FindCentersSpheres(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, const Vector3& local_center, vector<Vector3>& centers_tetra);
size_t MakeEnergy(const vector<Type>& Illum, const vector<Type>& squares, const Type scuare_surface, vector<Type>& energy);
Type IntegarteDirection(const int num_cell, const vector<Type>& Illum, const vector<Type>& directions);
size_t WriteFileSolution(const std::string NameFileOut, const std::vector<Type>& vector_energy, const vtkSmartPointer<vtkUnstructuredGrid>& UGrid);
size_t CenterOfTetra(const int NumberCell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Vector3& PointInTetra);
int FromSphericalToDecart(const int number_cur, const vector<Type>& all_directions, Vector3& direction);
int NormalsToCell(vtkCell* cell_tetra, Eigen::Matrix<Type, 4, 3>& normals);
int MakeFileDirectionsCenterTriangle(const std::string name_file_sphere_direction, const vtkSmartPointer<vtkUnstructuredGrid>& spehere_directions);
size_t ReadSphereDirectionDecartToSpherical(const std::string name_file_sphere_direction, vector<Type>& directions_all, vector<Type>& squares, Type& square_surface);
size_t ReadSphereDirectionDecartToSpherical(const std::string name_file_sphere_direction, vector<Type>& directions_all);
size_t TransformNetgenToVtkSurface(const std::string name_file_netgen, const std::string name_new_file_vtk);
Type IntegarteDirection(const vector<Type>& squares, const Type scuare_surface);
int FromPlaneToTetra(const Eigen::Matrix3d& inverse_transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& plane_coord, Eigen::Vector3d& tetra_coord);

template<typename T, typename T1>
T Min(const T a, const T1 b) {
	if (a < b) return a;
	return b;
}
#endif