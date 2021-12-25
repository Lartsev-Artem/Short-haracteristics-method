#include "short_characteristics_build_graph.h"

const double R = 0.5;
const Vector3 center(0, 0, 0);

int NewInitBoundarySet(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::set<IntId>& boundary_cells,
	std::set<IntId>& inner_boundary_faces);
int NewInitFacesState(const std::vector<IntId>& all_pairs_id, std::vector<IntId>& faces_state, const std::set<IntId>& inner_boundary_faces);
IntId NewGetNextCell(const Vector3& direction, const std::set<IntId>& boundary_id, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces);

int OldMainBuildGraphs(int argc, char* argv[]) {

	std::string main_file_direction = "D:\\Desktop\\FilesCourse\\";

	std::string name_file_settings = main_file_direction + "FOOsettings_file.txt";
	std::string name_file_grid = main_file_direction + "Sphere2519.vtk";//"Sphere109.vtk";//"MyVTK1.vtk";//
	std::string name_file_sphere_direction = main_file_direction + "SphereDir660.txt";
	std::string name_file_graph = main_file_direction + "TestClaster\\graph";//"GraphDirections\\graph";

	ReadStartSettings(name_file_settings, main_file_direction, name_file_sphere_direction, name_file_graph);

	std::vector<Type> directions;
	ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions);

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	ReadGridVtk(name_file_grid, unstructured_grid);
	std::cout << "Grid has cells: " << unstructured_grid->GetNumberOfCells() << '\n';


	std::vector<IntId> all_pairs_id;
	FindNeighborsPairFace(unstructured_grid, all_pairs_id);


	std::vector<IntId> faces_state;
	std::set<IntId> set_boundary_cells;
	std::set<IntId> set_graph;

	Vector3 direction;

	for (size_t cur_direction = 0; cur_direction < directions.size()/2; cur_direction++)
	{
		set_graph.clear();
		InitBoundarySet(unstructured_grid, set_boundary_cells);
		InitFacesState(all_pairs_id, faces_state);
		FromSphericalToDecart(cur_direction, directions, direction);

		//direction = Vector3(1, 0, 0);
	
		std::ofstream ofile;
		ofile.open(name_file_graph + to_string(cur_direction) + ".txt");
		if (!ofile.is_open()) {
			std::cout << "Error open file\n";
			std::cout << "file_graph is not opened for writing\n";
			return 1;
		}

		//std::vector<int> buf;
		while (set_boundary_cells.size()) {

			IntId id_cell;// = GetNextCell(direction, set_boundary_cells, unstructured_grid, faces_state, all_pairs_id);

			set_graph.emplace(id_cell);
			//buf.push_back(id_cell);
			ofile << id_cell << '\n';

			ReBuildSetBondary(id_cell, direction, set_boundary_cells, unstructured_grid, faces_state, all_pairs_id, set_graph);
		}
		ofile.close();

		std::cout << "Graph construction in the direction #" << cur_direction << " is completed\n";
		//WriteFileBoundary(main_file_direction + "Sphere565boundVTK.vtk", name_file_graph + "0.txt", unstructured_grid);
		//return -666;
	}
			
	return 0;
}

int MainBuildGraphs(int argc, char* argv[]) {

	std::string main_file_direction = "D:\\Desktop\\FilesCourse\\";

	std::string name_file_settings = main_file_direction + "FOOsettings_file.txt";
	std::string name_file_grid = main_file_direction + "MySphere.vtk";//"Sphere109.vtk";//"MyVTK1.vtk";//
	std::string name_file_sphere_direction = main_file_direction + "SphereDir660.txt";
	std::string name_file_graph = main_file_direction + "TestClaster\\graph";//"GraphDirections\\graph";

	ReadStartSettings(name_file_settings, main_file_direction, name_file_sphere_direction, name_file_graph);

	std::vector<Type> directions;
	ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions);

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	ReadGridVtk(name_file_grid, unstructured_grid);
	std::cout << "Grid has cells: " << unstructured_grid->GetNumberOfCells() << '\n';


	std::vector<IntId> all_pairs_id;
	FindNeighborsPairFace(unstructured_grid, all_pairs_id);


	std::vector<IntId> faces_state;
	std::set<IntId> set_boundary_cells;
	std::set<IntId> set_graph;

	std::set<IntId> set_inner_boundary_faces;

	Vector3 direction;

	for (size_t cur_direction = 0; cur_direction < directions.size() / 2; cur_direction++)
	{
		set_graph.clear();
		NewInitBoundarySet(unstructured_grid, set_boundary_cells, set_inner_boundary_faces);
		NewInitFacesState(all_pairs_id, faces_state, set_inner_boundary_faces);
		FromSphericalToDecart(cur_direction, directions, direction);

		direction = Vector3(1, 0, 0);

		std::ofstream ofile;
		ofile.open(name_file_graph + to_string(cur_direction) + ".txt");
		if (!ofile.is_open()) {
			std::cout << "Error open file\n";
			std::cout << "file_graph is not opened for writing\n";
			return 1;
		}

		//std::vector<int> buf;
		while (set_boundary_cells.size()) {

			IntId id_cell = NewGetNextCell(direction, set_boundary_cells, unstructured_grid, faces_state, all_pairs_id,set_inner_boundary_faces);

			static int n = 0;
			if (id_cell == -1) {
				if (n++ < 10) continue;
				else {
					ofile.close();
					std::cout << "Error.Complete " << set_graph.size() << " cells\n";
					break;
				}
			}
			n = 0;
			set_graph.emplace(id_cell);
			//buf.push_back(id_cell);
			ofile << id_cell << '\n';

			ReBuildSetBondary(id_cell, direction, set_boundary_cells, unstructured_grid, faces_state, all_pairs_id, set_graph);
		}
		ofile.close();

		std::cout << "Graph construction in the direction #" << cur_direction << " is completed\n";
		WriteFileBoundary(main_file_direction + "MySphereGraph.vtk", name_file_graph + "0.txt", unstructured_grid);
		return -666;
	}

	return 0;
}

int FindIdCellInBoundary(const Vector3& direction, const std::set<IntId>& inner_bound, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	const int cur_cell, const int cur_face) {

	// вершины треугольника.
	Vector3 P1;
	Vector3 P2;
	Vector3 P3;
	u_grid->GetCell(cur_cell)->GetFace(cur_face % 4)->GetPoints()->GetPoint(0, P1.data());
	u_grid->GetCell(cur_cell)->GetFace(cur_face % 4)->GetPoints()->GetPoint(1, P2.data());
	u_grid->GetCell(cur_cell)->GetFace(cur_face % 4)->GetPoints()->GetPoint(2, P3.data());



	Vector3 vertex;
	{
		// пока одна точка (центр)

		Type a = (P1 - P2).norm();
		Type b = (P3 - P2).norm();
		Type c = (P1 - P3).norm();
		Type S = a + b + c;

		vertex[0] = (a * P3[0] + b * P1[0] + c * P2[0]) / S;
		vertex[1] = (a * P3[1] + b * P1[1] + c * P2[1]) / S;
		vertex[2] = (a * P3[2] + b * P1[2] + c * P2[2]) / S;
	}

	// ищем пересечения vectrex->direction c гранями внутренней границы
	Vector3 intersect_point;
	int count = 0;
	int result = -1;
	for (auto& in_id: inner_bound) {	
		IntersectionWithPlane(u_grid->GetCell(in_id / 4)->GetFace(in_id % 4), vertex, direction, intersect_point);
		if (InTriangle(in_id / 4, u_grid, u_grid->GetCell(in_id / 4), in_id % 4, intersect_point))
			if(in_id!=cur_face){
				count++;
				result = in_id;
				//break;
			}
	}

	if (count > 1)
		std::cout << "Inner intersections: " << count << " ??\n";
	if(!count)
		std::cout << "No intersections: " << "\n";

	return result;
}
bool NewCheckCell(const IntId id_cell, const Vector3& direction, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces) {


	IntId face_state[4];
	//FindInAndOutFaces(direction, u_grid->GetCell(id_cell), face_state);
	FindInAndOutFaces(direction, id_cell, u_grid, face_state);

	for (size_t i = 0; i < 4; i++)
		if (face_state[i]) {// входящая грань

			int id_face = id_cell * 4 + i;
			if (faces_state[id_face] == 0 && all_pairs_id[id_face] == -1) {//  неопределенная внутренняя граница
				if (inner_boundary_faces.count(id_face)!=0) {
					int try_id = FindIdCellInBoundary(direction, inner_boundary_faces, u_grid, id_cell, id_face);
					if (try_id == -1) {
						std::cout << "NewCheckCell: wtf??\n";
						continue;
					}

					if (faces_state[try_id] == 1) { // если грань на другом конце определена
						faces_state[id_face] = 1;
						continue;
					}
					else
						return false;
				}
			}

			if (faces_state[id_cell * 4 + i] == 0) {// && faces_state[all_pairs_id[id_cell * 4 + i]] == 0) {
				// грани не определены
				return false;
			}
		}
	return true;
}

int FindIdCellInBoundary2(const Vector3& direction, const std::set<IntId>& inner_bound, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	const int cur_cell, const int cur_face, int* id) {

	// вершины треугольника.
	Vector3 P1;
	Vector3 P2;
	Vector3 P3;
	u_grid->GetCell(cur_cell)->GetFace(cur_face % 4)->GetPoints()->GetPoint(0, P1.data());
	u_grid->GetCell(cur_cell)->GetFace(cur_face % 4)->GetPoints()->GetPoint(1, P2.data());
	u_grid->GetCell(cur_cell)->GetFace(cur_face % 4)->GetPoints()->GetPoint(2, P3.data());




	// середины сторон (противолежащий точке P1,P2,P3 соответственно)
	Vector3 P11= (P3 + P2) / 2;
	Vector3 P22= (P3 + P1) / 2;
	Vector3 P33= (P2 + P1) / 2;

	// точки на медианах
	Vector3 vertex1 = P1 + (P11 - P1) / 3;
	Vector3 vertex2 = P2 + (P22 - P2) / 3;
	Vector3 vertex3 = P3 + (P33 - P3) / 3;


	// ищем пересечения vectrex->direction c гранями внутренней границы
	Vector3 intersect_point;
	int count = 0;

	for (auto& in_id : inner_bound) {
		IntersectionWithPlane(u_grid->GetCell(in_id / 4)->GetFace(in_id % 4), vertex1, direction, intersect_point);
		if (InTriangle(in_id / 4, u_grid, u_grid->GetCell(in_id / 4), in_id % 4, intersect_point))
			if (in_id != cur_face) {
				count++;
				id[0] = in_id;
				break;
			}
	}

	for (auto& in_id : inner_bound) {
		IntersectionWithPlane(u_grid->GetCell(in_id / 4)->GetFace(in_id % 4), vertex2, direction, intersect_point);
		if (InTriangle(in_id / 4, u_grid, u_grid->GetCell(in_id / 4), in_id % 4, intersect_point))
			if (in_id != cur_face) {
				count++;
				id[1] = in_id;
				break;
			}
	}

	for (auto& in_id : inner_bound) {
		IntersectionWithPlane(u_grid->GetCell(in_id / 4)->GetFace(in_id % 4), vertex3, direction, intersect_point);
		if (InTriangle(in_id / 4, u_grid, u_grid->GetCell(in_id / 4), in_id % 4, intersect_point))
			if (in_id != cur_face) {
				count++;
				id[2] = in_id;
				break;
			}
	}



	return 0;
}
bool NewCheckCell2(const IntId id_cell, const Vector3& direction, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces) {


	IntId face_state[4];
	//FindInAndOutFaces(direction, u_grid->GetCell(id_cell), face_state);
	FindInAndOutFaces(direction, id_cell, u_grid, face_state);

	for (size_t i = 0; i < 4; i++)
		if (face_state[i]) {// входящая грань

			int id_face = id_cell * 4 + i;
			if (faces_state[id_face] == 0 && all_pairs_id[id_face] == -1) {//  неопределенная внутренняя граница
				if (inner_boundary_faces.count(id_face) != 0) {
					int try_id[3] = { -1,-1,-1};
					FindIdCellInBoundary2(direction, inner_boundary_faces, u_grid, id_cell, id_face, try_id);
					if (try_id[0] == -1 || try_id[1] == -1 || try_id[2] == -1 ) {
						std::cout << "NewCheckCell: wtf??\n";
						continue;
					}

					if (faces_state[try_id[0]] == 1 && faces_state[try_id[1]] == 1 && faces_state[try_id[2]] == 1) { // если грань на другом конце определена
						faces_state[id_face] = 1;
						continue;
					}
					else
						return false;
				}
			}

			if (faces_state[id_cell * 4 + i] == 0) {// && faces_state[all_pairs_id[id_cell * 4 + i]] == 0) {
				// грани не определены
				return false;
			}
		}
	return true;
}
bool CheckCell(const IntId id_cell, const Vector3& direction, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id) {


	IntId face_state[4];
	//FindInAndOutFaces(direction, u_grid->GetCell(id_cell), face_state);
	FindInAndOutFaces(direction, id_cell, u_grid, face_state);

	for (size_t i = 0; i < 4; i++)
		if (face_state[i]) {
			if (faces_state[id_cell * 4 + i] == 0){// && faces_state[all_pairs_id[id_cell * 4 + i]] == 0) {
				// грани не определены
				return false;
			}
		}
	return true;
}

int InitBoundarySet(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::set<IntId>& boundary_cells) {

	boundary_cells.clear();
	int N = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkIdList> idc = vtkSmartPointer<vtkIdList>::New();

	for (size_t i = 0; i < N; ++i) {

		for (size_t j = 0; j < 4; ++j) {
			unstructured_grid->GetCellNeighbors(i, unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds(), idc);
			if (idc->GetNumberOfIds() == 0) {
				boundary_cells.emplace(i);
				break;
			}
		}
	}

	return 0;
}
int InitFacesState(const std::vector<IntId>& all_pairs_id, std::vector<IntId>& faces_state) {
	
	faces_state.assign(all_pairs_id.size(), 0);
	for (size_t i = 0; i < all_pairs_id.size(); i++) {
		if (all_pairs_id[i] == -1)
			faces_state[i] = 1;
	}
	return 0;
}
int ReadGridVtk(const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {
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
	return 0;
}
int ReBuildSetBondary(const IntId id_cell, const Vector3& direction, std::set<IntId>& boundary_id, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id, std::set<IntId>& set_graph) {

	IntId face_state[4];
	FindInAndOutFaces(direction, id_cell, u_grid, face_state);
	//FindInAndOutFaces(direction, u_grid->GetCell(id_cell), face_state);

	for (size_t i = 0; i < 4; i++)
		if (face_state[i] == 0) {
			int num_face = all_pairs_id[id_cell * 4 + i];

			int buf_s = set_graph.size();
			set_graph.emplace(num_face / 4);
			if (set_graph.size() != buf_s) {
				boundary_id.emplace(num_face / 4);
				set_graph.erase(num_face / 4);
			}
		}
	boundary_id.erase(id_cell);

	return 0;
}
size_t WriteFileBoundary(const std::string name_file_out, const std::string name_file_graph, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid) {

	int n = u_grid->GetNumberOfCells();

	vtkSmartPointer<vtkIntArray> bound_array =
		vtkSmartPointer<vtkIntArray>::New();

	bound_array->SetNumberOfTuples(n);

	ifstream ifile;
	ifile.open(name_file_graph);
	if (!ifile.is_open()) {
		std::cout << "Error open file\n";
		std::cout << "file_graph is not opened for reading\n";
		return 1;
	}

	int i = 0, el;
	while (ifile >> el)
		bound_array->SetTuple1(el, i++);
	ifile.close();

	vtkSmartPointer<vtkUnstructuredGrid> ungrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	ungrid = u_grid;
	ungrid->GetCellData()->SetActiveScalars("energy");
	ungrid->GetCellData()->SetScalars(bound_array);


	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileName(name_file_out.c_str());
	writer->SetInputData(ungrid);
	writer->Write();
	return 0;
}


int NewInitBoundarySet(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::set<IntId>& boundary_cells,
	std::set<IntId>& inner_boundary_faces) {

	boundary_cells.clear();
	inner_boundary_faces.clear();

	int N = unstructured_grid->GetNumberOfCells();

	vtkSmartPointer<vtkIdList> idc = vtkSmartPointer<vtkIdList>::New();

	for (size_t i = 0; i < N; ++i) {

		for (size_t j = 0; j < 4; ++j) {
			unstructured_grid->GetCellNeighbors(i, unstructured_grid->GetCell(i)->GetFace(j)->GetPointIds(), idc);
			if (idc->GetNumberOfIds() == 0) {

				Vector3 P(unstructured_grid->GetCell(i)->GetFace(j)->GetPoints()->GetPoint(0));
				if (P.norm()>R) { // внешняя сфера
					boundary_cells.emplace(i);
				}
				else {
					inner_boundary_faces.emplace(i * 4 + j);
					boundary_cells.emplace(i);
				}
				break;
			}
		}
	}

	return 0;
}
int NewInitFacesState(const std::vector<IntId>& all_pairs_id, std::vector<IntId>& faces_state, const std::set<IntId>& inner_boundary_faces) {

	faces_state.assign(all_pairs_id.size(), 0);
	for (size_t i = 0; i < all_pairs_id.size(); i++) {
		if (all_pairs_id[i] == -1)
			faces_state[i] = 1;
	}

	for (auto& el : inner_boundary_faces)
		faces_state[el] = 0;
	
	return 0;
}

IntId NewGetNextCell(const Vector3& direction, const std::set<IntId>& boundary_id, const vtkSmartPointer<vtkUnstructuredGrid>& u_grid,
	std::vector<IntId>& faces_state, const std::vector<IntId>& all_pairs_id, const std::set<IntId>& inner_boundary_faces) {
	// возвращает (случайный)индекс ячейки со всеми определенными входящими гранями

	for (auto el : boundary_id) {
		if (NewCheckCell(el, direction, u_grid, faces_state, all_pairs_id, inner_boundary_faces)) {
			for (size_t i = 0; i < 4; i++) {
				faces_state[el * 4 + i] = 1;
				if (all_pairs_id[el * 4 + i] != -1)
					faces_state[all_pairs_id[el * 4 + i]] = 1;
			}
			return el;
		}
	}

	std::cout << "NextCell is -1\n";
	return -1;
}