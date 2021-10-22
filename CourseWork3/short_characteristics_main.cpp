#include "short_characteristics_main.h"

Vector3 start_point_plane_coord; // начало координат плоскости
Matrix3 transform_matrix;   // матрица перехода из базового тетраэдра в плоскость
Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

Matrix3	straight_face;  // 3 узла интерпол€ции
Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра


int TestGetNodes(const int num_cur_out_face, const Eigen::Matrix4d& vertex_tetra) {

	Vector3 x;

	Vector3 node;

	switch (num_cur_out_face)
	{
	case 0:
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = 0;
			node[1] = straight_face.row(num_node)[0];
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выход€щей грани
			cout << "Number: " << num_cur_out_face << "\n";
			cout << x << '\n';

		}
		break;
	case 1:
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = 0;
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выход€щей грани
			cout << "Number: " << num_cur_out_face << "\n";
			cout << x << '\n';
		}
		break;
	case 2:
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = straight_face.row(num_node)[1];
			node[2] = 0;
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			cout << "Number: " << num_cur_out_face << "\n";
			cout << x << '\n';
		}// x->координата узла на выход€щей грани		}

		break;
	case 3:
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = inclined_face.row(num_node)[0];
			node[1] = inclined_face.row(num_node)[1];
			node[2] = 0;
			FromPlaneToTetra(inverse_transform_matrix, start_point_plane_coord, node, node);
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			
			cout << "Number: " << num_cur_out_face << "\n";
			cout << x << '\n';

		}
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

	
	//MakeFileDirectionsCenterTriangle(main_file_direction + "SphereDir.txt", unstructured_grid);

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
	SetVertexMatrix(0, unstructured_grid, vertex_tetra);
	cout << "vertex  " << vertex_tetra << '\n';
	for(int i=0;i<4;i++)
	TestGetNodes(i, vertex_tetra);
	return 0;

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
	for (int num_direction = 0; num_direction < count_direction; ++num_direction) {

		FromSphericalToDecart(num_direction, directions, direction);

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

			// перенос точки с выход€щей грани на вход€щую(получение точки на вход€щей грани)
			{			

				for (size_t num_out_face = 0; num_out_face < 4; ++num_out_face) {
					if (!face_state[num_out_face])  // выход€щие грани
						GetNodesValues(num_cell, unstructured_grid->GetCell(num_cell), num_out_face, face_state, direction, vertex_tetra, nodes_value, all_pairs_face,
							density, absorp_coef, rad_en_loose_rate,
							straight_face, inclined_face, transform_matrix, inverse_transform_matrix, start_point_plane_coord);
				}

				Type I_k_dir = GetValueInCenterCell(num_cell, unstructured_grid->GetCell(num_cell), centers_tetra[num_cell], direction, vertex_tetra, nodes_value, all_pairs_face,
					density, absorp_coef, rad_en_loose_rate,
					straight_face, inclined_face, transform_matrix, start_point_plane_coord);

				Illum[num_direction * count_direction + num_cell] = I_k_dir;
			}

		}
		/*---------------------------------- конец FOR по €чейкам----------------------------------*/

	}
	/*---------------------------------- конец FOR по направлени€м----------------------------------*/

	vector<Type> energy(unstructured_grid->GetNumberOfCells());
	MakeEnergy(Illum, squares, square_surface, energy);

	WriteFileSolution(name_file_vtk, energy, unstructured_grid);

	return 0;
}

