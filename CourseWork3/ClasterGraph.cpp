#include <stdio.h>
#include <iomanip>

#include <math.h>

#include<fstream>
#include<iostream>
#include<vector>
#include<map>
#include<set>
#include<string>

#include<omp.h>

#define FIGHT_GRID
using namespace std;
typedef double Type;
typedef int IntId ;
#define PI 3.14159265358979323846

#include<vtk-9.0/vtkUnstructuredGrid.h>

#ifdef FIGHT_GRID
#include <map>
#include <eigen3/Eigen/Dense>

Type r_inner = 0.5;
Eigen::Vector3d center_point(0, 0, 0);
struct Cell {
	Type P[3];
	Type P1[3];
	Type P2[3];
};

std::map<int, Cell> inner_bound;  //*global id face, vertex of face*/
std::set<int> inner_out_id;

int ReadInnerBoundary(const std::string name_file_inner_boundary, std::map<int, Cell>& inner_bound_loc) {

	ifstream ifile;

	ifile.open(name_file_inner_boundary);
	if (!ifile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}

	int N;
	ifile >> N;

	int id;
	Cell cell;
	while (ifile >> id >>
		cell.P[0] >> cell.P[1] >> cell.P[2] >>
		cell.P1[0] >> cell.P1[1] >> cell.P1[2] >>
		cell.P2[0] >> cell.P2[1] >> cell.P2[2]) {

		inner_bound_loc.emplace(id, cell);
	}


	ifile.close();
	return 0;
}
#endif // FIGHT_GRID

#ifdef vtkUnstructuredGrid_h
#include <vtk-9.0\vtkGenericDataObjectReader.h>
#include <vtk-9.0\vtkIdList.h>
#include <vtk-9.0\vtkCellData.h>
#include <vtk-9.0\vtkGenericDataObjectWriter.h>



int ReadPairsId(const std::string name_file_direction, vector<int>& all_pairs_id);

size_t NormalAndSquareFace(size_t NumberCell, size_t NumberFace, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Eigen::Vector3d& n) {
	vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetFace(NumberFace)->GetPointIds();

	Type P0[3], P1[3], P2[3];
	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);
	unstructuredgrid->GetPoint(idp->GetId(2), P2);

	Type a[3], b[3];
	for (size_t i = 0; i < 3; i++) {
		a[i] = P1[i] - P0[i];
		b[i] = P2[i] - P0[i];
	}
	n[0] = a[1] * b[2] - a[2] * b[1];
	n[1] = -a[0] * b[2] + a[2] * b[0];
	n[2] = a[0] * b[1] - a[1] * b[0];

	n.normalize();

	vtkSmartPointer<vtkIdList> idp2 = unstructuredgrid->GetCell(NumberCell)->GetPointIds();

	size_t id;
	for (size_t i = 0; i < 4; i++) {
		int count = 0;
		for (size_t j = 0; j < 3; j++)
			if (idp2->GetId(i) != idp->GetId(j))
				count++;
		if (count == 3) {
			id = i;
			break;
		}

	}

	Type sum = 0;
	Type P3[3];
	unstructuredgrid->GetPoint(idp2->GetId(id), P3);
	/*for (size_t i = 0; i < 3; i++){
		sum += n[i] * (P3[i] - P0[i]);
	}*/

	sum = P1[0] * (P2[1] - P3[1]) * P0[2] + P0[0] * (P3[1] - P2[1]) * P1[2] +
		P0[0] * (P1[1] - P3[1]) * P2[2] + P2[2] * (P1[0] * P3[1] - P1[0] * P0[1]) +
		P3[0] * (P0[2] * (P1[1] - P2[1]) + P1[2] * (P2[1] - P0[1]) + P2[2] * (P0[1] - P1[1]))
		+ P3[2] * (P1[0] * (P0[1] - P2[1]) + P0[0] * (P2[1] - P1[1])) +
		P2[0] * (P0[2] * (P3[1] - P1[1]) + P1[2] * (P0[1] - P3[1]) + P3[2] * (P1[1] - P0[1]));

	if (sum < 0)
		for (size_t i = 0; i < 3; i++)
			n[i] *= -1;
	return 0;
}

int WriteNormalFile(const std::string name_file_normals, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {
	ofstream ofile;

	ofile.open(name_file_normals);
	if (!ofile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}

	int n = unstructured_grid->GetNumberOfCells();
	Eigen::Vector3d normal;
	ofile << n << '\n';
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			NormalAndSquareFace(i, j, unstructured_grid, normal);
			ofile << normal[0] << ' ' << normal[1] << ' ' << normal[2] << '\n';
		}
		
	}
	ofile.close();
	return 0;
}

int TransformFromVtkToMyGrid(const std::string name_file_vtk, const std::string name_file_grid, 
	vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
		vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader_vtk->ReadAllScalarsOn();
	reader_vtk->SetFileName(name_file_vtk.c_str());
	reader_vtk->Update();

	if (reader_vtk->IsFileUnstructuredGrid()) {
		unstructured_grid = reader_vtk->GetUnstructuredGridOutput();
		unstructured_grid->Modified();
	}
	else {
		std::cout << "Error read file\n";
		std::cout << "file_vtk is not UnstructuredGrid\n";
		return 1;
	}

	ofstream ofile;

	ofile.open(name_file_grid);
	if (!ofile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}

	 int n = unstructured_grid->GetNumberOfPoints();

	Type P[3];
	for (size_t i = 0; i < n; i++)
		unstructured_grid->GetPoint(i, P);
		ofile << P[0] << ' ' << P[1] << ' ' << P[2] << '\n';

		n = unstructured_grid->GetNumberOfCells();
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < 3; j++)
				ofile << unstructured_grid->GetCell(i)->GetPointIds()->GetId(j) << ' ';
			ofile << unstructured_grid->GetCell(i)->GetPointIds()->GetId(3) << '\n';
		}

		//n *= 4;
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < 4; j++) {
				ofile << unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds()->GetId(0) << ' ';
				ofile << unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds()->GetId(1) << ' ';
				ofile << unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds()->GetId(2) << '\n';
			}
		}

	ofile.close();
	return 0;
	return 0;
}

int WritePairsId(const std::string name_file_pairs, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {
	
	auto GetNumberNeighborFace{[] (const int a, const int b, const int c, vtkCell* neighbor_cell)
		{

			vtkIdList* idc;

			int x, y, z;
			for (int i = 0; i < 4; i++)
			{
				idc = neighbor_cell->GetFace(i)->GetPointIds();
				x = idc->GetId(0);
				y = idc->GetId(1);
				z = idc->GetId(2);

				if (a == x && b == y && c == z) return i;
				else if (a == x && b == z && c == y) return i;
				else if (a == y && b == x && c == z) return i;
				else if (a == y && b == z && c == x) return i;
				else if (a == z && b == x && c == y) return i;
				else if (a == z && b == y && c == x) return i;

			}
			return -2;
		} };
	
	int count_unique_face = 0;
	const int N = unstructured_grid->GetNumberOfCells();
	std::vector<int>all_pairs_face(N * 4);
	for (int i = 0; i < N * 4; i++)
		all_pairs_face[i] = -2;

	vtkSmartPointer<vtkIdList> idp = vtkSmartPointer< vtkIdList>::New();
	vtkSmartPointer<vtkIdList> idc = vtkSmartPointer< vtkIdList>::New();

	int id_a, id_b, id_c;
	for (vtkIdType num_cell = 0; num_cell < N; ++num_cell) {

		for (int num_face = 0; num_face < 4; ++num_face) {
			if (all_pairs_face[num_cell * 4 + num_face] != -2) continue;
			++count_unique_face;

			idp = unstructured_grid->GetCell(num_cell)->GetFace(num_face)->GetPointIds();
			id_a = idp->GetId(0);
			id_b = idp->GetId(1);
			id_c = idp->GetId(2);

			/*Может быть проблема с указателями на списки!*/
			unstructured_grid->GetCellNeighbors(num_cell, idp, idc);

			if (idc->GetNumberOfIds() == 1) {
				int id_neighbor_cell = idc->GetId(0);
				int id_neighbor_face = GetNumberNeighborFace(id_a, id_b, id_c, unstructured_grid->GetCell(id_neighbor_cell));
				all_pairs_face[num_cell * 4 + num_face] = id_neighbor_cell * 4 + id_neighbor_face;
				all_pairs_face[id_neighbor_cell * 4 + id_neighbor_face] = num_cell * 4 + num_face;
			}
			else if (idc->GetNumberOfIds() == 0)
				all_pairs_face[num_cell * 4 + num_face] = -1; // граничная ячейка
			else
				std::cout << "More than 1 neighbor????\n";
		}

	}


	ofstream ofile;

	ofile.open(name_file_pairs);
	if (!ofile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}

	ofile << all_pairs_face.size() << '\n';
	for (size_t i = 0; i < all_pairs_face.size(); i++)
	{
		ofile << all_pairs_face[i] << '\n';
	}
	ofile.close();
	return 0;
}

int WriteInitBoundarySet(const std::string name_file_boundary, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	std::set<int> boundary_cells;
	int N = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkIdList> idc = vtkSmartPointer<vtkIdList>::New();

	for (size_t i = 0; i < N; ++i) {

		for (size_t j = 0; j < 4; ++j) {
			unstructured_grid->GetCellNeighbors(i, unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds(), idc);
			if (idc->GetNumberOfIds() == 0) {
#ifdef FIGHT_GRID
				Eigen::Vector3d P(unstructured_grid->GetCell(i)->GetFace(j)->GetPoints()->GetPoint(0));
				if ((P - Eigen::Vector3d(0, 0, 0)).norm() > r_inner) {
					boundary_cells.emplace(i);
					break;
				}
#else
				boundary_cells.emplace(i);
				break;
#endif
			}
		}
	}


	ofstream ofile;

	ofile.open(name_file_boundary);
	if (!ofile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}

	for (auto el : boundary_cells)
		ofile << el << '\n';

	ofile.close();
	return 0;
}

#ifdef FIGHT_GRID

int WriteBoundary(const std::string name_file_inner_boundary, std::vector<int>& all_pairs_id, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid) {

	std::vector<int> inner_boundary;

	for (size_t i = 0; i < all_pairs_id.size(); i++)
	{
		if (all_pairs_id[i] == -1) {
			Eigen::Vector3d P(unstructuredgrid->GetCell(i / 4)->GetFace(i % 4)->GetPoints()->GetPoint(0));
			if ((P - center_point).norm() < r_inner) {
				inner_boundary.push_back(i);
			}
		}
	}

	ofstream ofile;

	ofile.open(name_file_inner_boundary);
	if (!ofile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}

	Type* P;
	ofile << inner_boundary.size() << '\n';
	for (size_t i = 0; i < inner_boundary.size(); i++)
	{
		ofile << inner_boundary[i] << ' ';
		P = unstructuredgrid->GetCell(inner_boundary[i] / 4)->GetFace(inner_boundary[i] % 4)->GetPoints()->GetPoint(0);
		ofile << P[0] << ' ' << P[1] << ' ' << P[2] << ' ';

		P = unstructuredgrid->GetCell(inner_boundary[i] / 4)->GetFace(inner_boundary[i] % 4)->GetPoints()->GetPoint(1);
		ofile << P[0] << ' ' << P[1] << ' ' << P[2] << ' ';

		P = unstructuredgrid->GetCell(inner_boundary[i] / 4)->GetFace(inner_boundary[i] % 4)->GetPoints()->GetPoint(2);
		ofile << P[0] << ' ' << P[1] << ' ' << P[2] << '\n';
	}
	ofile.close();
	return 0;
}
#endif // FIGHT_GRID

int BuildSetForClaster(const std::string name_file_vtk, const std::string name_file_grid, const std::string name_file_pairs,
	const std::string name_file_boundary, const std::string name_file_normals, const std::string name_file_boundary_inner) {
	
	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
		vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader_vtk->ReadAllScalarsOn();
	reader_vtk->SetFileName(name_file_vtk.c_str());
	reader_vtk->Update();

	if (reader_vtk->IsFileUnstructuredGrid()) {
		unstructured_grid = reader_vtk->GetUnstructuredGridOutput();
		unstructured_grid->Modified();
	}
	else {
		std::cout << "Error read file\n";
		std::cout << "file_vtk is not UnstructuredGrid\n";
		return 1;
	}
//	TransformFromVtkToMyGrid(name_file_vtk, name_file_grid, unstructured_grid);
	std::cout << "Grid has Cell: " << unstructured_grid->GetNumberOfCells() << '\n';
	WritePairsId(name_file_pairs, unstructured_grid);
	WriteInitBoundarySet(name_file_boundary, unstructured_grid);

	WriteNormalFile(name_file_normals, unstructured_grid);

#ifdef FIGHT_GRID

	std::vector<int> all_pairs_id11;
	ReadPairsId(name_file_pairs, all_pairs_id11);
	WriteBoundary(name_file_boundary_inner, all_pairs_id11, unstructured_grid);

#endif // FIGHT_GRID

	cout << "Build end\n";
	return 0;
}

#endif // vtkUnstructuredGrid_h

#ifdef FIGHT_GRID
int FindIdCellInBoundary(const Type* direction, const std::map<int, Cell>& inner_bound, const int cur, int* id);
#endif // FIGHT_GRID

 size_t FromDecartToSphere(const Type* decart, Type& fi, Type& theta) {
	 Type x = decart[0];
	 Type y = decart[1];

	 //cout << x * x + y * y + decart[2] * decart[2] << '\n';

	 theta = atan(sqrt(x * x + y * y) / decart[2]) + PI / 2;

	 if (x <= 0)
		 fi = atan(y / x) + PI;
	 else if (x > 0 && y < 0)
		 fi = atan(y / x) + 2 * PI;
	 else if (x > 0 && y >= 0)
		 fi = atan(y / x);

	 return 0;
 }
size_t ReadSphereDirectionDecartToSpherical(const std::string name_file_sphere_direction, vector<Type>& directions_all) {

	ifstream ifile;

	ifile.open(name_file_sphere_direction);
	if (!ifile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}
	int N = 0;
	ifile >> N;
	directions_all.resize(2 * N);

	Type buf_s;
	int i = 0;
	Type P[3];
	Type theta;
	Type fi;
	for (int i = 0; i < N; i++) {
		ifile >> buf_s;
		ifile >> P[0] >> P[1] >> P[2];
		FromDecartToSphere(P, fi, theta);
		directions_all[i] = theta;
		directions_all[N + i] = fi;
	}
	ifile >> buf_s;
	ifile.close();
	return 0;
}
int ReadPairsId(const std::string name_file_direction, vector<int>& all_pairs_id) {
	ifstream ifile;

	ifile.open(name_file_direction);
	if (!ifile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}
	int N = 0;
	ifile >> N;
	all_pairs_id.resize(N);

	for (size_t i = 0; i < N; i++)
		ifile >> all_pairs_id[i];


	ifile.close();
	return 0;
}
int ReadNormals(const std::string name_file_normals, std::vector<Type>& normals) {
	ifstream ifile;

	ifile.open(name_file_normals);
	if (!ifile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}

	int n;
	ifile >> n;
	n *= 12;
	normals.resize(n);
	for (size_t i = 0; i < n; i++)
	{
		ifile >> normals[i];
	}
	ifile.close();
	return 0;
}
int ReadInitBoundarySet(const std::string name_file_direction, std::set<IntId>& set_graph) {
	ifstream ifile;

	ifile.open(name_file_direction);
	if (!ifile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}

	set_graph.clear();

	IntId buf;
	while (ifile >> buf)
	{
		set_graph.emplace(buf);
	}

	ifile.close();
	return 0;
}

int FromSphericalToDecart(const int number_cur, const vector<Type>& all_directions, Type* direction) {
	Type theta = all_directions[number_cur];
	Type fi = all_directions[all_directions.size() / 2 + number_cur];

	direction[0] = sin(theta) * cos(fi);
	direction[1] = sin(theta) * sin(fi);
	direction[2] = cos(theta);
	return 0;
}

int FindInAndOutFaces(const Type* direction, const int NumberCell, const std::vector<Type>& normals, int* face_state) {
	//face_state  -0=> выходящая грань,  1=> входящая  face_state.size=4!!!  

	Type normal[3];

	for (size_t i = 0; i < 4; ++i) {

		for (size_t j = 0; j < 3; j++)
		{
			normal[j] = normals[NumberCell * 12 + i * 3 + j];
		}
		//Eigen::Vector3d normall;
		//NormalAndSquareFace(NumberCell, i, unstructuredgrid, normal);
		const Type eps = 1e-6;

		Type scalar = 0;
		for (size_t i = 0; i < 3; i++)
			scalar += normal[i] * direction[i];

		if (scalar < -eps)
			face_state[i] = 1;
		else if (scalar > eps)
			face_state[i] = 0;
		else 
			face_state[i] = 1;
#ifdef vtkUnstructuredGrid_h
		//else
			//std::cout << "FindInAndOutFaces: error directions\n";
#endif // vtkUnstructuredGrid_h
	}

	return 0;
}


bool CheckCell(const IntId id_cell, const Type* direction, const std::vector<Type>& normals,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id) {


	IntId face_state[4];
	//FindInAndOutFaces(direction, u_grid->GetCell(id_cell), face_state);
	FindInAndOutFaces(direction, id_cell, normals, face_state);
#ifdef FIGHT_GRID

	
	int id[3]; // номера пересечений (global)

	for (size_t i = 0; i < 4; i++)
		if (face_state[i]) {
			
			int id_face = id_cell * 4 + i;
			if (all_pairs_id[id_face] == -1 && faces_state[id_face] == 0&&false) { // внутренняя  сфера граница 
				// Find Cell если все ячейки для 3х вершин грани определены, то continue
				// трасировка из вершин по направлению (НУЖНА СЕТКА)
				FindIdCellInBoundary(direction, inner_bound, id_face, id);

				//if (id[0] != -1 || id[1] != -1 || id[2] != -1)
				//	cout << "id\n";
				for (int k = 0; k < 1; k++) {   // ЗДЕСЬ НЕ 0-1 А 0-3
					if (id[k] != -1)
						if (faces_state[id[k]] == 0)
							return false;
				}
				//inner_out_id.erase(id_cell);
				return true;

			}
			if (faces_state[id_cell * 4 + i] == 0)
				return false;
		}
#else
	for (size_t i = 0; i < 4; i++)
		if (face_state[i]) {
			if (faces_state[id_cell * 4 + i] == 0) {// && faces_state[all_pairs_id[id_cell * 4 + i]] == 0) {
				// грани не определены
				return false;
			}
		}
#endif
	return true;
}
IntId GetNextCell(const Type* direction, const std::set<IntId>& boundary_id, const std::vector<Type>& normals,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id) {
	// возвращает (случайный)индекс ячейки со всеми определенными входящими гранями

	for (auto el : boundary_id) {
		if (CheckCell(el, direction, normals, faces_state, all_pairs_id)) {
			
			{
				for (size_t i = 0; i < 4; i++) {
					faces_state[el * 4 + i] = 1;
					if (all_pairs_id[el * 4 + i] != -1)
						faces_state[all_pairs_id[el * 4 + i]] = 1;
				}
			}
			return el;
		}
	}

	std::cout << "NextCell is -1\n";
	return -1;
}

int InitFacesState(const std::vector<IntId>& all_pairs_id, std::vector<IntId>& faces_state) {

	faces_state.assign(all_pairs_id.size(), 0);
	for (size_t i = 0; i < all_pairs_id.size(); i++) {
		if (all_pairs_id[i] == -1)
			faces_state[i] = 1;
	}


#ifdef FIGHT_GRID
	for (auto &el : inner_bound)
		faces_state[el.first] = 0;
#endif
	return 0;
}


int ReBuildSetBondary(const IntId id_cell, const Type* direction, std::set<IntId>& boundary_id, const std::vector<Type>& normals,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id, std::set<IntId>& set_graph) {

	IntId face_state[4];
	FindInAndOutFaces(direction, id_cell, normals, face_state);
	//FindInAndOutFaces(direction, u_grid->GetCell(id_cell), face_state);

#ifndef FIGHT_GRID
	for (size_t i = 0; i < 4; i++)
		if (face_state[i] == 0) {
			int num_face = all_pairs_id[id_cell * 4 + i];
			if (num_face != -1) {
				int buf_s = set_graph.size();
				set_graph.emplace(num_face / 4);
				if (set_graph.size() != buf_s) {
					boundary_id.emplace(num_face / 4);
					set_graph.erase(num_face / 4);
				}
			}
		}
	boundary_id.erase(id_cell);
#else
	for (size_t i = 0; i < 4; i++)
		if (face_state[i] == 0) {
			int num_face = all_pairs_id[id_cell * 4 + i];
			if (num_face != -1) {
				int buf_s = set_graph.size();
				set_graph.emplace(num_face / 4);
				if (set_graph.size() != buf_s) {
					boundary_id.emplace(num_face / 4);
					set_graph.erase(num_face / 4);
				}
			}
			else {
				if (inner_bound.count(id_cell * 4 + i))
				{
					//добавить элементы на другом конце сферы
					{
						int id[3];
						FindIdCellInBoundary(direction, inner_bound, id_cell * 4 + i, id);
						if (id[0] != -1) {
							//if (faces_state[id[0]] == 0)
								int buf_s = set_graph.size();
							set_graph.emplace(id[0] / 4);
							if (set_graph.size() != buf_s) {
								boundary_id.emplace(id[0] / 4);
								faces_state[id[0]] = 1;
								set_graph.erase(id[0] / 4);
								inner_out_id.erase(id[0] / 4);
							}
						}
					}
					
					//return 0;// не удалять ячейку с границы если она на врнутренней границе
				}
			}
		}
	//если это внутрення граница сразу не удалять!
	
	boundary_id.erase(id_cell);
#endif
	return 0;
}


size_t ReadStartSettings(const std::string name_file_settings,  std::string& name_file_sphere_direction, std::string& name_file_graph, std::string& name_file_boundary,
	std::string& name_file_grid_pairs, std::string&  name_file_grid_normals, std::string& name_file_errors) {

	std::ifstream ifile;
	ifile.open(name_file_settings);
	if (!ifile.is_open()) {
		std::cerr << " Error : file settings is not open !\n";
		return 1;
	}

	std::string str; // переменная для перевода строки при чтении из файла

	getline(ifile, name_file_sphere_direction);
	getline(ifile, name_file_graph);
	getline(ifile, name_file_boundary);
	getline(ifile, name_file_grid_pairs);
	getline(ifile, name_file_grid_normals);
	getline(ifile, name_file_errors);

	ifile.close();
	return 0;
}

int WriteStringToFile(const std::string name_file, char* str, std::ios::openmode mode) {
	
	std::ofstream ofile;
	ofile.open(name_file, mode); //file_settings
	if (!ofile.is_open()) {
		std::cout << "Error open file\n";
		std::cout << "WriteStringToFile is not opened for writing\n";
	}
	ofile << str;
	ofile.close();
	return 0;
}

int mainWork(int argc, char* argv[]) {

#ifndef vtkUnstructuredGrid_h
	std::string main_file_direction = "D:\\Desktop\\FilesCourse\\TestClaster\\";

	std::string name_file_grid = "D:\\Desktop\\FilesCourse\\MySphere.vtk";
	std::string name_file_sphere_direction = "D:\\Desktop\\FilesCourse\\SphereDir660.txt";
	std::string name_file_graph = main_file_direction + "graph";

	std::string  name_file_boundary = main_file_direction + "Bounadary.txt";
	std::string  name_file_grid_pairs = main_file_direction + "Pairs.txt";
	std::string  name_file_grid_normals = main_file_direction + "Normals.txt";
	std::string  name_file_boundary_inner = main_file_direction + "BounadaryInner.txt";



	BuildSetForClaster(name_file_grid, main_file_direction + "grid.txt", name_file_grid_pairs, name_file_boundary, name_file_grid_normals, name_file_boundary_inner);
	return 0;
#else
	const std::string  file_settings = "D:\\Desktop\\FilesCourse\\TestClaster\\settings_claster_file.txt";
	std::string  name_file_sphere_direction;
	std::string  name_file_graph;

	std::string  name_file_boundary;
	std::string  name_file_grid_pairs;
	std::string  name_file_grid_normals;
	//std::string  name_file_boundary_inner;
	std::string name_file_errors;
	ReadStartSettings(file_settings, name_file_sphere_direction, name_file_graph, name_file_boundary, name_file_grid_pairs, name_file_grid_normals,
		name_file_errors);
#endif

	Type full_time = -omp_get_wtime();

	std::vector<Type> directions;
	ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions);


	std::vector<int> all_pairs_id;
	ReadPairsId(name_file_grid_pairs, all_pairs_id);

	std::vector<Type> normals;
	ReadNormals(name_file_grid_normals, normals);

	omp_set_num_threads(4);
#pragma omp parallel default(none) shared(all_pairs_id,directions,normals)  
	{
		std::vector<IntId> faces_state;
		std::set<IntId> set_boundary_cells;
		std::set<IntId> set_graph;

		Type direction[3];

#pragma omp for 
		for (int cur_direction = 0; cur_direction < directions.size() / 2; cur_direction++) {

			set_graph.clear();
			ReadInitBoundarySet(name_file_boundary, set_boundary_cells);
			InitFacesState(all_pairs_id, faces_state);
			FromSphericalToDecart(cur_direction, directions, direction);

			std::ofstream ofile;
			ofile.open(name_file_graph + to_string(cur_direction) + ".txt");
			if (!ofile.is_open()) {
				std::cout << "Error open file\n";
				std::cout << "file_graph is not opened for writing\n";
				break;
				//return 1;
			}

			std::vector<int> buf;
			while (set_boundary_cells.size()) {

				IntId id_cell = GetNextCell(direction, set_boundary_cells, normals, faces_state, all_pairs_id);

				if (id_cell < 0) {
					ofile.close();

#ifndef vtkUnstructuredGrid_h
					char str[60];
					sprintf_s(str, "Errors direction: %d, error cell: %d \n", cur_direction, id_cell);
					WriteStringToFile(name_file_errors, str, std::ios::app);  //name_file_errors
#endif
					break;
					//return 1;
				}

				set_graph.emplace(id_cell);
				buf.push_back(id_cell); 
#pragma omp critical
				{
					ofile << id_cell << '\n';
					ofile << fixed;
				}

				ReBuildSetBondary(id_cell, direction, set_boundary_cells, normals, faces_state, all_pairs_id, set_graph);
			}
			ofile.close();
			printf("Graph construction in the direction # %d is completed\n", cur_direction);

#ifdef vtkUnstructuredGrid_h
			printf("Graph construction in the direction # %d is completed\n", cur_direction);
#endif

			//WriteFileBoundary(main_file_direction + "Sphere565boundVTK.vtk", name_file_graph + "0.txt", unstructured_grid);
			//return -666;
		}
	}

	full_time += omp_get_wtime();
	char str[30];
	sprintf_s(str, "\nFull time work: %lf\n", full_time);

#ifdef vtkUnstructuredGrid_h
	printf(str);
#else
	WriteStringToFile(file_settings, str, std::ios::app);  //file_settings
#endif
	return 0;
}


#ifdef EIGEN_MACROS_H
size_t SetBasis(const Type* start_point, Eigen::Vector3d& normal, Eigen::Matrix3d& basis) {
	/*по начальной точке и нормале строит локальный базис картинной плоскости (vec1, vec2).
	  нормаль дана. задаем один вектор произвольно(ортагонально normal). третий вектор из векторного произведения*/
	Eigen::Vector3d vec_1;
	Eigen::Vector3d vec_2;

	if (fabs(normal[1]) < 1e-20) {
		vec_1[0] = 0;
		vec_1[1] = 1;
		vec_1[2] = 0;
	}
	else {
		vec_1[0] = 1;
		vec_1[2] = 0;
		vec_1[1] = -(normal[0] * vec_1[0] + normal[2] * vec_1[2]) / normal[1];  //св-во скалярного произведения (N, vec1)==0
	}

	// правельная ориентация базиса плоскости
	if (normal[1] < 0)
		for (int i = 0; i < 3; ++i)
			vec_1[i] *= -1;

	// обычное векторное умножение. Eigen временно и не нужен!!!
	Eigen::Vector3d c = normal.cross(vec_1);

	for (size_t i = 0; i < 3; ++i)
		vec_2[i] = -c(i);

	vec_1.normalize();
	vec_2.normalize();

	basis.row(0) = vec_1;
	basis.row(1) = vec_2;
	basis.row(2) = normal;

	return 0;
}
size_t Make2dPoint(const Type* start, const Eigen::Matrix3d& local_basis, const Type* point, Eigen::Vector3d& new_point) {



	for (size_t i = 0; i < 3; i++)
		new_point[i] = 0;

	//перевод 3d точки в 2d (в локальном базисе {start, local_basis}) 
	for (size_t k = 0; k < 3; k++) {
		new_point[0] += (point[k] - start[k]) * local_basis(0, k);
		new_point[1] += (point[k] - start[k]) * local_basis(1, k);
	}
	return 0;
}

#ifdef FIGHT_GRID

bool InTriangle(const int num_cell, const std::map<int, Cell>& inner_bound, const Type* XX, const Type* dir) {
	/*face --- треугольник, X --- точка для проверки*/

	// вершины треугольника
	Type AA[3] = { inner_bound.find(num_cell)->second.P[0],inner_bound.find(num_cell)->second.P[1],inner_bound.find(num_cell)->second.P[2] };
	Type BB[3] = { inner_bound.find(num_cell)->second.P1[0],inner_bound.find(num_cell)->second.P1[1],inner_bound.find(num_cell)->second.P1[2] };
	Type CC[3] = { inner_bound.find(num_cell)->second.P2[0],inner_bound.find(num_cell)->second.P2[1],inner_bound.find(num_cell)->second.P2[2] };
	
	Eigen::Vector3d A, B, C, X;
	{
		Type Xx[3] = { XX[0], XX[1], XX[2] };
		Eigen::Vector3d n = { dir[0],dir[1],dir[2] };
		Eigen::Matrix3d basis;
		SetBasis(AA, n, basis);
		//cout << basis << '\n';
		Make2dPoint(AA, basis, AA, A);
		Make2dPoint(AA, basis, BB, B);
		Make2dPoint(AA, basis, CC, C);
		Make2dPoint(AA, basis, Xx, X);
	}

	// линейная алгебра
	Type r1 = (A[0] - X[0]) * (B[1] - A[1]) - (B[0] - A[0]) * (A[1] - X[1]);
	Type r2 = (B[0] - X[0]) * (C[1] - B[1]) - (C[0] - B[0]) * (B[1] - X[1]);
	Type r3 = (C[0] - X[0]) * (A[1] - C[1]) - (A[0] - C[0]) * (C[1] - X[1]);

	Type eps = 1e-6;
	if (r1 < -eps && r2 < -eps && r3 < -eps)
		return true;
	else if (r1 > eps && r2 > eps && r3 > eps)
		return true;
	else return false;
}
int FindIdCellInBoundary(const Type* direction,  const std::map<int, Cell>& inner_bound, const int cur, int* id) {
	
	id[0] = -1;
	id[1] = -1;
	id[2] = -1;

	if (!inner_bound.count(cur)) {
		printf("Inner bound is empty??\n");
		return 0;
	}
	Type vertex_1[3] = { inner_bound.find(cur)->second.P[0],inner_bound.find(cur)->second.P[1],inner_bound.find(cur)->second.P[2] };
	Type vertex_2[3] = { inner_bound.find(cur)->second.P1[0],inner_bound.find(cur)->second.P1[1],inner_bound.find(cur)->second.P1[2] };
	Type vertex_3[3] = { inner_bound.find(cur)->second.P2[0],inner_bound.find(cur)->second.P2[1],inner_bound.find(cur)->second.P2[2] };
	// вершины треугольника плохо. Надо внутренние точки
	
	Type vertex[3];
	{
	// пока одна точка (центр)
		
		Type a = sqrt((vertex_1[0] - vertex_2[0]) + (vertex_1[1] - vertex_2[1]) + (vertex_1[2] - vertex_2[2]));
		Type b= sqrt((vertex_3[0] - vertex_2[0]) + (vertex_3[1] - vertex_2[1]) + (vertex_3[2] - vertex_2[2]));
		Type c = sqrt((vertex_1[0] - vertex_3[0]) + (vertex_1[1] - vertex_3[1]) + (vertex_1[2] - vertex_3[2]));
		Type S = a + b + c;

		vertex[0] = (a * vertex_3[0] + b * vertex_1[0] + c * vertex_2[0]) / S;
		vertex[1] = (a * vertex_3[1] + b * vertex_1[1] + c * vertex_2[1]) / S;
		vertex[2] = (a * vertex_3[2] + b * vertex_1[2] + c * vertex_2[2]) / S;
	}


	for (auto &el : inner_bound) {
		if (InTriangle(el.first, inner_bound, vertex, direction) && el.first != cur) { // проверка пересечения луча(vertex+dir)
			id[0] = el.first;
			break;
		}
	}

	for (auto &el : inner_bound) {
		if (InTriangle(el.first, inner_bound, vertex, direction) && el.first != cur) { // проверка пересечения луча(vertex+dir)
			id[1] = el.first;
			break;
		}
	}

	for (auto &el : inner_bound) {
		if (InTriangle(el.first, inner_bound, vertex, direction) && el.first != cur) { // проверка пересечения луча(vertex+dir)
			id[2] = el.first;
			break;
		}
	}

	return 0;
}

int FindOutInnerCell(const std::map<int, Cell>& inner_bound, const Type* dir, const std::vector<Type>& normals, set<int>& inner_out_id) {
	int face[4];
	inner_out_id.clear();
	for (auto& el : inner_bound)
	{
		FindInAndOutFaces(dir, el.first/4, normals, face);

		for (int i = 0; i < 4; i++){
			if (face[i] == 1)// входящая
				if (el.first%4 == i)// граница сферы равна тек. грани
					inner_out_id.emplace(el.first / 4);

		}

	}
	return 0;
}
#endif // FIGHT_GRID

#ifdef  vtkUnstructuredGrid_h
size_t WriteFileBoundary(const std::string name_file_out, const std::string name_file_graph, const std::string name_file_grid) {

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	{
		vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
			vtkSmartPointer<vtkGenericDataObjectReader>::New();
		reader_vtk->ReadAllScalarsOn();
		reader_vtk->SetFileName(name_file_grid.c_str());
		reader_vtk->Update();

		if (reader_vtk->IsFileUnstructuredGrid()) {
			unstructured_grid = reader_vtk->GetUnstructuredGridOutput();
			unstructured_grid->Modified();
		}
		else {
			std::cout << "Error read file\n";
			std::cout << "file_vtk is not UnstructuredGrid\n";
			return 1;
		}
	}
	int n = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkIntArray> bound_array =
		vtkSmartPointer<vtkIntArray>::New();

	bound_array->SetNumberOfTuples(n);
	for (size_t i = 0; i < n; i++)
	{
		bound_array->SetTuple1(i, 0);
	}

	ifstream ifile;
	ifile.open(name_file_graph);
	if (!ifile.is_open()) {
		std::cout << "Error open file\n";
		std::cout << "file_graph is not opened for reading\n";
		return 1;
	}

	string str;
	int i = 1, el;
	while (ifile >> el) {
			bound_array->SetTuple1(el, i++);
	}
	//ifile >> el;
	/*while (!ifile.eof()) {
		if(ifile >> el)
		bound_array->SetTuple1(el, i++);
		getline(ifile, str);
	}*/
	ifile.close();

	vtkSmartPointer<vtkUnstructuredGrid> ungrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	ungrid = unstructured_grid;
	ungrid->GetCellData()->SetActiveScalars("energy");
	ungrid->GetCellData()->SetScalars(bound_array);


	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileName(name_file_out.c_str());
	writer->SetInputData(ungrid);
	writer->Write();
	return 0;
}
size_t WriteFileBoundary(const std::string name_file_out, const std::set<int>& id, const std::string name_file_grid) {

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	{
		vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
			vtkSmartPointer<vtkGenericDataObjectReader>::New();
		reader_vtk->ReadAllScalarsOn();
		reader_vtk->SetFileName(name_file_grid.c_str());
		reader_vtk->Update();

		if (reader_vtk->IsFileUnstructuredGrid()) {
			unstructured_grid = reader_vtk->GetUnstructuredGridOutput();
			unstructured_grid->Modified();
		}
		else {
			std::cout << "Error read file\n";
			std::cout << "file_vtk is not UnstructuredGrid\n";
			return 1;
		}
	}
	int n = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkIntArray> bound_array =
		vtkSmartPointer<vtkIntArray>::New();

	bound_array->SetNumberOfTuples(n);
	for (size_t i = 0; i < n; i++)
	{
		bound_array->SetTuple1(i, 0);
	}



	for(auto el :id){
		bound_array->SetTuple1(el, 1);
	}
	

	vtkSmartPointer<vtkUnstructuredGrid> ungrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	ungrid = unstructured_grid;
	ungrid->GetCellData()->SetActiveScalars("energy");
	ungrid->GetCellData()->SetScalars(bound_array);


	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileName(name_file_out.c_str());
	writer->SetInputData(ungrid);
	writer->Write();
	return 0;
}
size_t WriteInnerSphere(const std::string name_file_out) {
	
	
	
	/*Type AA[3] = { inner_bound.find(num_cell)->second.P[0],inner_bound.find(num_cell)->second.P[1],inner_bound.find(num_cell)->second.P[2] };
	Type BB[3] = { inner_bound.find(num_cell)->second.P1[0],inner_bound.find(num_cell)->second.P1[1],inner_bound.find(num_cell)->second.P1[2] };
	Type CC[3] = { inner_bound.find(num_cell)->second.P2[0],inner_bound.find(num_cell)->second.P2[1],inner_bound.find(num_cell)->second.P2[2] };*/
	return 0;
}
#endif //  vtkUnstructuredGrid_h

int main(int argc, char* argv[]) {

	std::string main_file_direction = "D:\\Desktop\\FilesCourse\\TestClaster\\";

	std::string name_file_grid = "D:\\Desktop\\FilesCourse\\MySphere.vtk";
	std::string name_file_sphere_direction = "D:\\Desktop\\FilesCourse\\SphereDir660.txt";
	std::string name_file_graph = main_file_direction + "graphTest";

	std::string  name_file_boundary = main_file_direction + "Bounadary.txt";
	std::string  name_file_grid_pairs = main_file_direction + "Pairs.txt";
	std::string  name_file_grid_normals = main_file_direction + "Normals.txt";

	std::string  name_file_boundary_inner = main_file_direction + "BounadaryInner.txt"; 

	BuildSetForClaster(name_file_grid, main_file_direction + "grid.txt", name_file_grid_pairs, name_file_boundary, name_file_grid_normals, name_file_boundary_inner);

	std::vector<Type> directions;
	ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions);

	std::vector<int> all_pairs_id;
	ReadPairsId(name_file_grid_pairs, all_pairs_id);

#ifdef FIGHT_GRID
	ReadInnerBoundary(name_file_boundary_inner, inner_bound);
#endif // FIGHT_GRID


	std::vector<Type> normals;
	ReadNormals(name_file_grid_normals, normals);

	omp_set_num_threads(1);
#pragma omp parallel default(none) shared(all_pairs_id,directions,normals)  
	{
		std::vector<IntId> faces_state;
		std::set<IntId> set_boundary_cells;
		std::set<IntId> set_graph;

		Type direction[3];

#pragma omp for 
		for (int cur_direction = 0; cur_direction<1/*cur_direction < directions.size() / 2*/; cur_direction++)
		{
			set_graph.clear();
			ReadInitBoundarySet(name_file_boundary, set_boundary_cells);
			InitFacesState(all_pairs_id, faces_state);
			FromSphericalToDecart(cur_direction, directions, direction);
			direction[0] = 1;
			direction[1] = 0;
			direction[2] = 0;
			FindOutInnerCell(inner_bound, direction, normals, inner_out_id);
			//WriteFileBoundary(main_file_direction + "WTF1.vtk", inner_out_id, name_file_grid);

			std::ofstream ofile;
			ofile.open(name_file_graph + to_string(cur_direction) + ".txt");
			if (!ofile.is_open()) {
				std::cout << "Error open file\n";
				std::cout << "file_graph is not opened for writing\n";
				break;
				//return 1;
			}

			std::vector<int> buf;
			while (set_boundary_cells.size()) {

				IntId id_cell = GetNextCell(direction, set_boundary_cells, normals, faces_state, all_pairs_id);

				if (id_cell < 0) {
					ofile.close();
					break;
					//return 1;
				}

				set_graph.emplace(id_cell);
				buf.push_back(id_cell);
#pragma omp critical
				{
					ofile << id_cell << '\n';
					ofile << fixed;
				}

				ReBuildSetBondary(id_cell, direction, set_boundary_cells, normals, faces_state, all_pairs_id, set_graph);
			}
			ofile.close();

			printf("Graph construction in the direction # %d is completed\n", cur_direction);
			//std::cout << "Graph construction in the direction #" << cur_direction << " is completed\n";
			//WriteFileBoundary(main_file_direction + "Sphere565boundVTK.vtk", name_file_graph + "0.txt", unstructured_grid);
			//return -666;

		}
	}
	
	WriteFileBoundary(main_file_direction + "WTF.vtk",  name_file_graph + to_string(0) + ".txt", name_file_grid);
	return 0;
}

#endif