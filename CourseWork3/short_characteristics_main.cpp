#include "short_characteristics_main.h"

Vector3 start_point_plane_coord; // начало координат плоскости
Matrix3 transform_matrix;   // матрица перехода из базового тетраэдра в плоскость
Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

Matrix3	straight_face;  // 3 узла интерпол€ции
Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра

int PrintTransformTetra(const Eigen::Matrix4d& vertex_tetra, vtkCell* cell) {

	Vector3 PP;
	for (size_t i = 0; i < 4; i++)
	{
		Vector3 P(cell->GetPoints()->GetPoint(i));
		cout << "Global P" << i << "  " << P[0] << ' ' << P[1] << ' ' << P[2] << '\n';
		FromGlobalToLocalTetra(vertex_tetra, P, PP);
		cout << "Local P" << i << "  " << PP[0] << ' ' << PP[1] << ' ' << PP[2] << "\n\n";
	}

	return 0;
}
int PrintCurCell(vtkCell* cell) {

	static int num_cell = 0;
	Type P[3];

	cout << "Cell: " << num_cell++ << '\n';
	for (size_t i = 0; i < 4; i++)
	{
		cout << "Face: " << i << '\n';
		for (size_t j = 0; j < 3; j++)
		{
			cell->GetFace(i)->GetPoints()->GetPoint(j, P);
			cout << P[0] << ' ' << P[1] << ' ' << P[2] << '\n';
		}
		cout << '\n';

	}

	return 0;
}
size_t WriteFileSolution(const std::string NameFileOut, const std::vector<int>& sort_id, const vtkSmartPointer<vtkUnstructuredGrid>& UGrid);
size_t WriteFileSolutionIllum(const std::string NameFileOut, const std::vector<Type>& Illum, const vtkSmartPointer<vtkUnstructuredGrid>& UGrid);
int TestCalculateNodeValue(const int num_cur_cell, vtkCell* cur_cell, const int num_cur_out_face, const int* face_state,
	const Vector3& direction,
	std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face,
	const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x, const Eigen::Matrix<Type,4,3> normals) {

	Vector3 x0;

	for (size_t num_in_face = 0; num_in_face < 4; ++num_in_face) {
		if (!face_state[num_in_face]) continue;  // обрабатываем только входные грани

		
       if (direction.dot(normals.row(num_in_face)) < eps) {
			std::cout << "direction.normal==0, SOS!!!\n"; continue; // плоскость параллельна направлению
		} 

		if (nodes_value.find(Min(all_pairs_face[num_cur_cell * 4 + num_in_face], num_cur_cell * 4 + num_in_face))->second[0] < 0) {
			cout<< Min(all_pairs_face[num_cur_cell * 4 + num_in_face], num_cur_cell * 4 + num_in_face) << "  Undefine face!!!\n";
		}

		IntersectionWithPlane(cur_cell->GetFace(num_in_face), x, direction, x0);

		if (InTriangle(cur_cell->GetFace(num_in_face), x0)) {

			//cout << all_pairs_face[num_cur_cell * 4 + num_in_face] << "  " << num_cur_cell * 4 + num_in_face << '\n';
			int global_num_in_face = Min(all_pairs_face[num_cur_cell * 4 + num_in_face], num_cur_cell * 4 + num_in_face);
			Type s = (x - x0).norm();

			// значение на вход€щей грани
			Type I_x0 = CalculateIllumeOnInnerFace(num_in_face, global_num_in_face, vertex_tetra, x, x0, nodes_value,
				straight_face, inclined_face, transform_matrix, start_point_plane_coord);

			int global_num_out_face = Min(all_pairs_face[num_cur_cell * 4 + num_cur_out_face], num_cur_cell * 4 + num_cur_out_face);
			if (global_num_out_face == -1)global_num_out_face = num_cur_cell * 4 + num_cur_out_face;

			nodes_value.find(global_num_out_face)->second[num_node] = 10. + (I_x0 - 10.) * exp(-s);
			break;
		}

	}//for num_in_face

	return 0;
}
int TestGetNodes(const int num_cur_cell, vtkCell* cur_cell, const int num_cur_out_face, const Eigen::Matrix4d& vertex_tetra,
	const int* face_state,
	const Vector3& direction,
	std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face, const Eigen::Matrix<Type, 4, 3> normals) {

	Vector3 x;

	Vector3 node;

	switch (num_cur_out_face)
	{
		//1->0
	case 0:
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = 0;
			node[1] = straight_face.row(num_node)[0];
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выход€щей грани

			TestCalculateNodeValue(num_cur_cell, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				num_node, vertex_tetra, x, normals);
		}
		cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 1:
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = 0;
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выход€щей грани				
			TestCalculateNodeValue(num_cur_cell, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				num_node, vertex_tetra, x, normals);
		}
		cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 2:
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = straight_face.row(num_node)[1];
			node[2] = 0;
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			TestCalculateNodeValue(num_cur_cell, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				num_node, vertex_tetra, x, normals);
		}// x->координата узла на выход€щей грани		}
		cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 3:
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = inclined_face.row(num_node)[0];
			node[1] = inclined_face.row(num_node)[1];
			node[2] = 0;
			FromPlaneToTetra(inverse_transform_matrix, start_point_plane_coord, node, node);
			FromLocalToGlobalTetra(vertex_tetra, node, x);						
			TestCalculateNodeValue(num_cur_cell, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				num_node, vertex_tetra, x, normals);
		}
		cout << "Number: " << num_cur_out_face << "\n";
		break;
	default:
		std::cout << "Number face is not {0,1,2,3}????\n";
		break;
	}

	return 0;
}


int main(int argc, char* argv[])
{
	std::string name_file_settings;
	name_file_settings = "D:\\Desktop\\FilesCourse\\settings_file.txt";

	std::string main_file_direction = "D:\\Desktop\\FilesCourse\\";

	size_t class_file_vtk;
	std::string name_file_vtk;
	std::string name_file_sphere_direction;
	std::string out_file_grid_vtk;
	

	//TransformFileDecartToSphere(name_file_sphere_direction, name_file_sphere_direction+"1.txt");
	//TransformNetgenToVtk(main_file_direction + "Sphere109.txt", main_file_direction + "Sphere109.vtk");
	//return 0;
	//TransformNetgenToVtkSurface(main_file_direction + "SphereDirectionNetral.txt", main_file_direction + "SphereDir.vtk");
	
	if (ReadStartSettings(name_file_settings, class_file_vtk, name_file_vtk, name_file_sphere_direction, out_file_grid_vtk)) {

		std::cout << "Error reading the start settings\n";
		return 1;
	}

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	// скал€рные данные сетки (unstructured_grid)
	vtkDataArray* density;
	vtkDataArray* absorp_coef;
	vtkDataArray* rad_en_loose_rate;

	Type _clock = -omp_get_wtime();
	if (ReadFileVtk(class_file_vtk, name_file_vtk, unstructured_grid, density, absorp_coef, rad_en_loose_rate, true)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the vtk_grid file: " << _clock << "\n";

	
//	MakeFileDirectionsCenterTriangle(main_file_direction + "SphereDir.txt", unstructured_grid);

	// theta - fi
	vector<Type> directions;
	vector<Type> squares;
	Type square_surface;
	_clock = -omp_get_wtime();
	//ReadSphereDirectionVtk(1, name_file_sphere_direction, directions);
	ReadSphereDirectionDecartToSpherical(name_file_sphere_direction, directions, squares, square_surface);
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the sphere_direction file: " << _clock << "\n";

	//MakeFileDirections(main_file_direction + "Direction.txt", directions);
	//MakeFileCellData(main_file_direction + "CellData.txt", unstructured_grid);

	// ”пор€доченные индексы €чеек по данному направлению
	vector<int> sorted_id_cell(unstructured_grid->GetNumberOfCells());

	vector<Vector3> centers_tetra;

		_clock = -omp_get_wtime();
		//ReadCentersFromGeneralFile(main_file_direction + "CellData.txt", centers_tetra);
		FindCentersSpheres(unstructured_grid, centers_tetra);
		//ReadCenters(main_file_direction + "centers.txt", centers_tetra);
		_clock += omp_get_wtime();
		std::cout << "\n Reading time of the CellData file: " << _clock << "\n";
		

	InitGlobalValue(start_point_plane_coord, transform_matrix, inverse_transform_matrix, straight_face, inclined_face);


	Eigen::Matrix4d vertex_tetra;
	/* x1 x2 x3 x4
	*  y1 y2 y3 y4
	*  z1 z2 z3 z4
	*  1  1  1  1
	*/


	std::vector<int> all_pairs_face;
	_clock = -omp_get_wtime();
	int count_unique_face = FindNeighborsPairFace(unstructured_grid, all_pairs_face);  //сделать файл!!!!!!
	_clock += omp_get_wtime();
	std::cout << "\n Finding time of the all_pairs_face: " << _clock << "\n";


	std::map<int, Vector3> nodes_value;
	InitNodesValue(all_pairs_face, nodes_value);  // выделение пам€ти под узловые значени€(без инициализации)


	const int count_direction = directions.size() / 2;

	// значени€ излучени€ в центре €чеек (dir^k*N_dir+N_cell) \
	номер строки --- направление, столбец --- €чейка

	std::vector<Type> Illum(unstructured_grid->GetNumberOfCells() * count_direction);

	Vector3 direction; // главное направление

	/*---------------------------------- далее FOR по направлени€м----------------------------------*/
	//for (int num_direction = 0; num_direction < count_direction; ++num_direction) 
	{
		int num_direction = 0;
		FromSphericalToDecart(num_direction, directions, direction);

		direction = Vector3(1, 0, 0);

		_clock = -omp_get_wtime();
		SortCellsGrid(direction, unstructured_grid, sorted_id_cell, centers_tetra);
		_clock += omp_get_wtime();
		std::cout << "\n Sort" << num_direction << " time : " << _clock << '\n';


		Vector3 x;
		Vector3 x0;
		Eigen::Matrix<Type, 4, 3> normals;

		int num_cell;
		/*---------------------------------- далее FOR по €чейкам----------------------------------*/
		for (int h = 0; h < unstructured_grid->GetNumberOfCells(); ++h) {
			num_cell = sorted_id_cell[h];

			SetVertexMatrix(num_cell, unstructured_grid, vertex_tetra);

			int face_state[4];  // 0=> выход€ща€ грань,  1=> вход€ща€   
			FindInAndOutFaces(direction, unstructured_grid->GetCell(num_cell), face_state);

			NormalsToCell(unstructured_grid->GetCell(num_cell), normals);

			PrintCurCell(unstructured_grid->GetCell(num_cell));
			PrintTransformTetra(vertex_tetra, unstructured_grid->GetCell(num_cell));
			// перенос точки с выход€щей грани на вход€щую(получение точки на вход€щей грани)
			{			

				for (size_t num_out_face = 0; num_out_face < 4; ++num_out_face) {
					if (!face_state[num_out_face])  // выход€щие грани
						TestGetNodes(num_cell, unstructured_grid->GetCell(num_cell), num_out_face, vertex_tetra, face_state, direction, nodes_value, all_pairs_face,
							normals);
						//GetNodesValues(num_cell, unstructured_grid->GetCell(num_cell), num_out_face, face_state, direction, vertex_tetra, nodes_value, all_pairs_face, \
							density, absorp_coef, rad_en_loose_rate, \
							straight_face, inclined_face, transform_matrix, inverse_transform_matrix, start_point_plane_coord);
				}

				Type I_k_dir = GetValueInCenterCell(num_cell, unstructured_grid->GetCell(num_cell), centers_tetra[num_cell], direction, vertex_tetra, nodes_value, all_pairs_face,
					density, absorp_coef, rad_en_loose_rate,
					straight_face, inclined_face, transform_matrix, start_point_plane_coord);

				Illum[num_direction * count_direction + num_cell] = I_k_dir;
			}

		}
		/*---------------------------------- конец FOR по €чейкам----------------------------------*/

		std::cout << "End direction number: " << num_direction << '\n';
	}
	/*---------------------------------- конец FOR по направлени€м----------------------------------*/

	vector<Type> energy(unstructured_grid->GetNumberOfCells());
	MakeEnergy(Illum, squares, square_surface, energy);
	WriteFileSolutionIllum(out_file_grid_vtk, Illum, unstructured_grid);
	//WriteFileSolution(out_file_grid_vtk, energy, unstructured_grid);
//	WriteFileSolution(out_file_grid_vtk, sorted_id_cell, unstructured_grid);

	return 0;
}
size_t WriteFileSolutionIllum(const std::string NameFileOut, const std::vector<Type>& Illum, const vtkSmartPointer<vtkUnstructuredGrid>& UGrid) {

	int n = UGrid->GetNumberOfCells();

	vtkSmartPointer<vtkDoubleArray> IllumArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	IllumArray->SetNumberOfTuples(n);
	for (size_t i = 0; i < n; i++)
		IllumArray->SetTuple1(i, Illum[i]);

	vtkSmartPointer<vtkUnstructuredGrid> ungrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	ungrid = UGrid;
	ungrid->GetCellData()->SetActiveScalars("energy");
	ungrid->GetCellData()->SetScalars(IllumArray);


	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileName(NameFileOut.c_str());
	writer->SetInputData(ungrid);
	writer->Write();
	return 0;
}
