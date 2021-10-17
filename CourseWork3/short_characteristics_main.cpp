#include "short_characteristics_main.h"

Vector3 start_point_plane_coord; // начало координат плоскости
Matrix3 transform_matrix;   // матрица перехода из базового тетраэдра в плоскость
Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

Matrix3	straight_face;  // 3 узла интерпол€ции
Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

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
	//TransformNetgenToVtk(main_file_direction + "Sphere.txt", main_file_direction + "Sphere.vtk");

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
	std::cout << "\n Reading time of the vtk file: " << ((Type)_clock) / CLOCKS_PER_SEC << "\n";

	
	MakeFileCellData(main_file_direction + "CellData.txt", unstructured_grid);
	
	// theta - fi
	vector<Type> directions;
	ReadSphereDirectionVtk(2, name_file_sphere_direction, directions);

	std::cout << "Integrate: " << IntegarteDirection(directions)/(4*PI) << '\n';

	vector<int> sorted_id_cell(unstructured_grid->GetNumberOfCells());
	vector<Vector3> centers_tetra(unstructured_grid->GetNumberOfCells());
	ReadCentersFromGeneralFile(main_file_direction + "CellData.txt", centers_tetra);

	Vector3 direction(1, 1, 0); // главное направление (тест)
	_clock = -omp_get_wtime();
	SortCellsGrid(direction, unstructured_grid, sorted_id_cell, centers_tetra);
	_clock += omp_get_wtime();
	std::cout << "\n Sort time : " << ((Type)_clock) / CLOCKS_PER_SEC << "\n";


	std::vector<int> all_pairs_face;
	int count_unique_face = FindNeighborsPairFace(unstructured_grid, all_pairs_face);

	std::map<int, Vector3> nodes_value;
	SetNodesValue(all_pairs_face, nodes_value);

	// значени€ излучени€ в центре €чеек (dir^k*N_dir+N_cell)
	std::vector<Type> Illum(unstructured_grid->GetNumberOfCells() * directions.size());


	// 3 узла интерпол€ции
	{
		straight_face << 1. / 6, 1. / 6, 1,
			2. / 3, 1. / 6, 1,
			1. / 6, 2. / 3, 1;
	}
	
	// 3 узла интерпол€ции на наклонной плоскости
	{ 
		inclined_face <<
			0, sqrt(2. / 3), 1,
			sqrt(2) / 4, 1. / (2 * sqrt(6)), 1,
			-sqrt(2) / 4, 1. / (2 * sqrt(6)), 1;
	}

	//ћатрицы перехода из стандартного тетраэдра в координаты наклонной плоскости 
	{ transform_matrix <<
		-1. / sqrt(2), 1. / sqrt(2), 0,
		-1. / sqrt(6), -1. / sqrt(6), sqrt(2. / 3),
		1. / sqrt(3), 1. / sqrt(3), 1. / sqrt(3);
	}

	//ћатрицы перехода из наклонной плоскости в  координаты стандартного тетраэдра
	{
		inverse_transform_matrix <<
			-1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
			1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
			0, sqrt(2. / 3), 1. / sqrt(3);
	}

	// Ќачало координата плоскости
	start_point_plane_coord << 0.5, 0.5, 0;
	

	Eigen::Matrix4d vertex_tetra; 
	/* x1 x2 x3 x4
	*  y1 y2 y3 y4
	*  z1 z2 z3 z4
	*  1  1  1  1
	*/
	
	int num_cell=0;

    SetVertexMatrix(num_cell, unstructured_grid, vertex_tetra);
	
	int face_state[4];  // 0=> выход€ща€ грань,  1=> вход€ща€   
	FindInAndOutFaces(direction, unstructured_grid->GetCell(num_cell), face_state);


	// перенос точки с выход€щей гране на вход€щую(получение точки на вход€щей грани)
	{
		Vector3 x;
		Vector3 x0;
		Eigen::Matrix<Type, 4, 3> normals;
		Vector3 center;

		ifstream ifile; ifile.open(main_file_direction + "CellData.txt");
		ReadNextCellData(ifile, normals, center);
			ifile.close();

		for (size_t num_out_face = 0; num_out_face < 4; ++num_out_face){
			if (!face_state[num_out_face])  // выход€щие грани
				
				GetNodesValues(num_cell, unstructured_grid->GetCell(num_cell), num_out_face, face_state, direction, vertex_tetra, nodes_value, all_pairs_face,
					density, absorp_coef, rad_en_loose_rate,
					straight_face, inclined_face, transform_matrix, start_point_plane_coord);
		}

		Type I_k_dir = GetValueInCenterCell(num_cell, unstructured_grid->GetCell(num_cell), center, direction, vertex_tetra, nodes_value, all_pairs_face,
			density, absorp_coef, rad_en_loose_rate, 
			straight_face, inclined_face, transform_matrix, start_point_plane_coord);

	}

	return 0;
}

