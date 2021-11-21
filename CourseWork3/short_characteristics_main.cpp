#include "short_characteristics_main.h"
typedef int IntId;

Vector3 start_point_plane_coord; // начало координат плоскости
Matrix3 transform_matrix;   // матрица перехода из базового тетраэдра в плоскость
Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

Matrix3	straight_face;  // 3 узла интерпол€ции
Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

Vector3 center_local_sphere;  // центр описанной сферы около стандартного тетраэдра

int GetAverageSizeCell(const std::string name_file_vtk) {

	auto SquareFace{ [](const size_t NumberCell, const size_t NumberFace, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid) {
		
		vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetFace(NumberFace)->GetPointIds();

		
		Type n1,n2,n3;
		Type P0[3], P1[3], P2[3];
		unstructuredgrid->GetPoint(idp->GetId(0), P0);
		unstructuredgrid->GetPoint(idp->GetId(1), P1);
		unstructuredgrid->GetPoint(idp->GetId(2), P2);

		Type a[3], b[3];
		for (size_t i = 0; i < 3; i++) {
			a[i] = P1[i] - P0[i];
			b[i] = P2[i] - P0[i];
		}
		
		n1 = a[1] * b[2] - a[2] * b[1];
		n2 = -a[0] * b[2] + a[2] * b[0];
		n3 = a[0] * b[1] - a[1] * b[0];

		return sqrt(pow(n1, 2) + pow(n2, 2) + pow(n3, 2)) / 2;

	} };

	auto GetVolumeCell{ [](const size_t NumberCell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid) {
		Type V = 0;
		Type P0[3], P1[3], P2[3], P3[3];

		vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetPointIds();
		unstructuredgrid->GetPoint(idp->GetId(0), P0);
		unstructuredgrid->GetPoint(idp->GetId(1), P1);
		unstructuredgrid->GetPoint(idp->GetId(2), P2);
		unstructuredgrid->GetPoint(idp->GetId(3), P3);

		Type a[3], b[3], c[3];
		for (size_t i = 0; i < 3; i++) {
			a[i] = P1[i] - P0[i];
			b[i] = P2[i] - P0[i];
			c[i] = P3[i] - P0[i];
		}

		V = a[0] * (b[1] * c[2] - c[1] * b[2]) - a[1] * (b[0] * c[2] - c[0] * b[2]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
		return fabs(V) / 6;
	} };


	vtkSmartPointer<vtkUnstructuredGrid> unstructuredgrid = 
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	{
		vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
			vtkSmartPointer<vtkGenericDataObjectReader>::New();
		reader_vtk->ReadAllScalarsOn();
		reader_vtk->SetFileName(name_file_vtk.c_str());
		reader_vtk->Update();

		if (reader_vtk->IsFileUnstructuredGrid()) {
			unstructuredgrid = reader_vtk->GetUnstructuredGridOutput();
			unstructuredgrid->Modified();
		}
		else {
			std::cout << "Error read file\n";
			std::cout << "file_vtk is not UnstructuredGrid\n";
			return 1;
		}
	}

	Type averageH = 0;
	Type minH = 1000;
	const size_t NumberOfCells = unstructuredgrid->GetNumberOfCells();

	Type VolumeCell = 0;
	Type S0 = 0;
	Type S1 = 0;
	Type S2 = 0;
	Type S3 = 0;

	for (size_t curCell = 0; curCell < NumberOfCells; ++curCell) {

		// 4  площади к текущей €чейке
		S0 = SquareFace(curCell, 0, unstructuredgrid);
		S1 = SquareFace(curCell, 1, unstructuredgrid);
		S2 = SquareFace(curCell, 2, unstructuredgrid);
		S3 = SquareFace(curCell, 3, unstructuredgrid);

		// объем текущей €чейки		
		VolumeCell = GetVolumeCell(curCell, unstructuredgrid);

		//–адиус вписанной сферы
		Type curH = 3 * VolumeCell / (S0 + S1 + S2 + S3);

		averageH += curH;
		if (curH < minH)
			minH = curH;
	}

	std::cout << fixed << setprecision(16) << "averageH= " << averageH / NumberOfCells << '\n';
	std::cout << fixed << setprecision(16) << "MinH= " << minH << '\n';
	return 0;
}

int GetFunction(const std::string name_file_vtk, const std::string name_file_out);
int ReadGraph(const std::string name_file_graph, std::vector<IntId>& sorted_id_cell) {
	ifstream ifile;
	ifile.open(name_file_graph);
	if (!ifile.is_open()) {
		std::cout << "Error open file\n";
		std::cout << "file_graph is not opened for reading\n";
		return 1;
	}
	
	IntId i;
	int count = 0;
	while (ifile >> i) {
		sorted_id_cell[count++] = i;
	}
	ifile.close();
	return 0;
}

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
int TestCalculateNodeValue(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& UGrid,vtkCell* cur_cell, 
	const int num_cur_out_face, const int* face_state,
	const Vector3& direction,
	std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face,
	const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x, const Eigen::Matrix<Type,4,3> normals) {

	Vector3 x0;

	for (size_t num_in_face = 0; num_in_face < 4; ++num_in_face) {
		if (!face_state[num_in_face]) continue;  // обрабатываем только входные грани

		Type buf = fabs(direction.dot(normals.row(num_in_face)));
       if (fabs(direction.dot(normals.row(num_in_face))) < eps) {
			std::cout << "direction.normal==0, SOS!!!\n"; continue; // плоскость параллельна направлению
		} 

		IntersectionWithPlane(cur_cell->GetFace(num_in_face), x, direction, x0);

		if (InTriangle(num_cur_cell,UGrid,cur_cell, num_in_face, x0)) {

			//cout << all_pairs_face[num_cur_cell * 4 + num_in_face] << "  " << num_cur_cell * 4 + num_in_face << '\n';
			int global_num_in_face = Min(all_pairs_face[num_cur_cell * 4 + num_in_face], num_cur_cell * 4 + num_in_face);
			
			Type s = (x - x0).norm();

			// значение на вход€щей грани
			Type I_x0 = CalculateIllumeOnInnerFace(num_in_face, global_num_in_face, vertex_tetra, x, x0, nodes_value,
				straight_face, inclined_face, transform_matrix, start_point_plane_coord);

			//cout <<"X0 VALUE " << I_x0 << '\n';

			int global_num_out_face = Min(all_pairs_face[num_cur_cell * 4 + num_cur_out_face], num_cur_cell * 4 + num_cur_out_face);
			if (global_num_out_face == -1)global_num_out_face = num_cur_cell * 4 + num_cur_out_face;


		//	cout << "X VALUE " << 10. + (I_x0 - 10.) * exp(-s) << '\n';
			Type Ie = 10;
			Type k = 10;
			if (x.norm() > 0.3) { Ie = 0; k = 1; }

			Type I;
			if (s>1e-10)
				 I = Ie * (1 - exp(-s*k)) + I_x0 * exp(-s*k);
			else
				 I = Ie * (1+s*k) - I_x0 * s*k;
				

			nodes_value.find(global_num_out_face)->second[num_node] = I;// Ie* (1 - exp(-s)) + I_x0 * exp(-s); //Ie + (I_x0 - Ie) * exp(-s);
			
			if (I < 0)
				nodes_value.find(global_num_out_face)->second[num_node] = 0;
			break;
		}

	}//for num_in_face

	return 0;
}
int TestGetNodes(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& UGrid, vtkCell* cur_cell, const int num_cur_out_face, const Eigen::Matrix4d& vertex_tetra,
	const int* face_state,
	const Vector3& direction,
	std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face, const Eigen::Matrix<Type, 4, 3> normals) {

	Vector3 x;

	Vector3 node;

	switch (num_cur_out_face)
	{
		
	case 1:// 1->2
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = 0;
			node[1] = straight_face.row(num_node)[0];
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выход€щей грани

			TestCalculateNodeValue(num_cur_cell, UGrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				num_node, vertex_tetra, x, normals);
		}
		//cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 2://2->0
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = 0;
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выход€щей грани				
			TestCalculateNodeValue(num_cur_cell, UGrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				num_node, vertex_tetra, x, normals);
		}
		//cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 0: //0->3
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = straight_face.row(num_node)[1];
			node[2] = 0;
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			TestCalculateNodeValue(num_cur_cell, UGrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				num_node, vertex_tetra, x, normals);
		}// x->координата узла на выход€щей грани		}
	//	cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 3: //3->1
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = inclined_face.row(num_node)[0];
			node[1] = inclined_face.row(num_node)[1];
			node[2] = 0;
			FromPlaneToTetra(inverse_transform_matrix, start_point_plane_coord, node, node);
			FromLocalToGlobalTetra(vertex_tetra, node, x);						
			TestCalculateNodeValue(num_cur_cell, UGrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				num_node, vertex_tetra, x, normals);
		}
	//	cout << "Number: " << num_cur_out_face << "\n";
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
	std::string name_file_graph = main_file_direction + "GraphDirections\\Direction660\\Sphere565\\graph"; //"TestClaster\\graph";//"GraphDirections\\graph";
	std::string out_file_grid_vtk;


	//TransformFileDecartToSphere(name_file_sphere_direction, name_file_sphere_direction+"1.txt");
	//TransformNetgenToVtk(main_file_direction + "Grids\\Sphere.txt", main_file_direction + "Sphere565.vtk");
	//return 0;
	
	
	if (ReadStartSettings(name_file_settings, class_file_vtk, name_file_vtk, name_file_sphere_direction, out_file_grid_vtk)) {

		std::cout << "Error reading the start settings\n";
		return 1;
	}


	//MainBuildGraphs(argc, argv);
	//return 0;

	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	// скал€рные данные сетки (unstructured_grid)
	vtkDataArray* density;
	vtkDataArray* absorp_coef;
	vtkDataArray* rad_en_loose_rate;

	//TransformNetgenToVtkSurface(main_file_direction + "SphereDirectionNetral.txt", main_file_direction + "SphereDir.vtk");
	//TransformNetgenToVtk(main_file_direction + "MySphere.txt", main_file_direction + "MySphere.vtk");
	//return 0;

	Type _clock = -omp_get_wtime();
	if (ReadFileVtk(class_file_vtk, name_file_vtk, unstructured_grid, density, absorp_coef, rad_en_loose_rate, true)) {
		std::cout << "Error reading the file vtk\n";
		return 1;
	}
	_clock += omp_get_wtime();
	std::cout << "\n Reading time of the vtk_grid file: " << _clock << "\n";


	//MakeFileDirectionsCenterTriangle(main_file_direction + "SphereDir.txt", unstructured_grid);

	//return 0;

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
	vector<IntId> sorted_id_cell(unstructured_grid->GetNumberOfCells());

	//vector<Vector3> centers_tetra;

		_clock = -omp_get_wtime();
		//ReadCentersFromGeneralFile(main_file_direction + "CellData.txt", centers_tetra);
		//FindCentersSpheres(unstructured_grid, centers_tetra);
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
	for (int num_direction = 0; num_direction < count_direction; ++num_direction) 
	{
	//	int num_direction = 0;
		FromSphericalToDecart(num_direction, directions, direction);

		//direction = Vector3(1, 0, 0);

		_clock = -omp_get_wtime();
		//SortCellsGrid(direction, unstructured_grid, sorted_id_cell, centers_tetra);
		ReadGraph(name_file_graph + to_string(num_direction) + ".txt", sorted_id_cell);
		_clock += omp_get_wtime();
		std::cout << "\n Sort" << num_direction << " time : " << _clock << '\n';

	//	std::reverse(sorted_id_cell.begin(), sorted_id_cell.end());

		Vector3 x;
		Vector3 x0;
		Eigen::Matrix<Type, 4, 3> normals;

		int num_cell;
		/*---------------------------------- далее FOR по €чейкам----------------------------------*/
		for (int h = 0; h < unstructured_grid->GetNumberOfCells(); ++h) {
			num_cell = sorted_id_cell[h];

			SetVertexMatrix(num_cell, unstructured_grid, vertex_tetra);

			int face_state[4];  // 0=> выход€ща€ грань,  1=> вход€ща€   
			//FindInAndOutFaces(direction, unstructured_grid->GetCell(num_cell), face_state);
			FindInAndOutFaces(direction, num_cell, unstructured_grid, face_state);

			NormalsToCell(num_cell, unstructured_grid, normals);
		//	NormalsToCell(unstructured_grid->GetCell(num_cell), normals);

			if ((num_cell == 2114 || num_cell == 2178) &&false) {
					PrintCurCell(unstructured_grid->GetCell(num_cell));
					PrintTransformTetra(vertex_tetra, unstructured_grid->GetCell(num_cell));
					// перенос точки с выход€щей грани на вход€щую(получение точки на вход€щей грани)
			}
			{			

				for (size_t num_out_face = 0; num_out_face < 4; ++num_out_face) {
					if (!face_state[num_out_face])  // выход€щие грани
						TestGetNodes(num_cell, unstructured_grid, unstructured_grid->GetCell(num_cell), num_out_face, vertex_tetra, face_state, direction, 
							nodes_value, all_pairs_face,
							normals);
						//GetNodesValues(num_cell, unstructured_grid->GetCell(num_cell), num_out_face, face_state, direction, vertex_tetra, nodes_value, all_pairs_face, \
							density, absorp_coef, rad_en_loose_rate, \
							straight_face, inclined_face, transform_matrix, inverse_transform_matrix, start_point_plane_coord);
				}
				Vector3 center;
				CenterOfTetra(num_cell, unstructured_grid, center);
				Type I_k_dir = GetValueInCenterCell(num_cell, unstructured_grid, unstructured_grid->GetCell(num_cell), center, direction, vertex_tetra, nodes_value, all_pairs_face,
					density, absorp_coef, rad_en_loose_rate,
					straight_face, inclined_face, transform_matrix, start_point_plane_coord);

				Illum[num_direction * unstructured_grid->GetNumberOfCells() + num_cell] = I_k_dir;

				/*int buf = num_cell;
				if (all_pairs_face[num_cell] < num_cell && all_pairs_face[num_cell] != -1) buf = all_pairs_face[num_cell];
				Illum[num_direction * unstructured_grid->GetNumberOfCells() + num_cell] = nodes_value.find(buf)->second[0];*/ //I_k_dir;
			}

		}
		/*---------------------------------- конец FOR по €чейкам----------------------------------*/

		std::cout << "End direction number: " << num_direction << '\n';
	}
	/*---------------------------------- конец FOR по направлени€м----------------------------------*/

	vector<Type> energy(unstructured_grid->GetNumberOfCells());
	MakeEnergy(Illum, squares, square_surface, energy);
	//WriteFileSolutionIllum(out_file_grid_vtk, Illum, unstructured_grid);
	WriteFileSolution(out_file_grid_vtk, energy, unstructured_grid);
//	WriteFileSolution(out_file_grid_vtk, sorted_id_cell, unstructured_grid);

	GetFunction(out_file_grid_vtk, main_file_direction + "E.txt");

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


bool InTraceTriangle(const int num_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cell_face, int number_face, const Eigen::Vector3d& XX) {
	/*face --- треугольник, X --- точка дл€ проверки*/

	// вершины треугольника
	Type AA[3];
	Type BB[3];
	Type CC[3];
	cell_face->GetFace(number_face)->GetPoints()->GetPoint(0, AA);
	cell_face->GetFace(number_face)->GetPoints()->GetPoint(1, BB);
	cell_face->GetFace(number_face)->GetPoints()->GetPoint(2, CC);


	Eigen::Vector3d A, B, C, X;
	{
		Type Xx[3] = { XX[0], XX[1], XX[2] };
		Eigen::Vector3d n = { 1./sqrt(3),1./sqrt(3),1. / sqrt(3) };
		Eigen::Matrix3d basis;
		SetBasis(Xx, n, basis);
		//cout << basis << '\n';
		Make2dPoint(Xx, basis, AA, A);
		Make2dPoint(Xx, basis, BB, B);
		Make2dPoint(Xx, basis, CC, C);
		Make2dPoint(Xx, basis, Xx, X);
	}

	// линейна€ алгебра
	Type r1 = (A[0] - X[0]) * (B[1] - A[1]) - (B[0] - A[0]) * (A[1] - X[1]);
	Type r2 = (B[0] - X[0]) * (C[1] - B[1]) - (C[0] - B[0]) * (B[1] - X[1]);
	Type r3 = (C[0] - X[0]) * (A[1] - C[1]) - (A[0] - C[0]) * (C[1] - X[1]);

	if (r1 < 0 && r2 < 0 && r3 < 0)
		return true;
	else if (r1 > 0 && r2 > 0 && r3 > 0)
		return true;
	else return false;
}
int GetIdTrace(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<int>& id_cells) {
	
	for (size_t i = 0; i < unstructured_grid->GetNumberOfCells(); i++)
	{
		for (size_t j = 0; j < 4; j++)
		{
			if (InTraceTriangle(i, unstructured_grid, unstructured_grid->GetCell(i), j, Eigen::Vector3d(-1, -1, -1))) {
				id_cells.push_back(i);
				break;
			}
		}
	}

	return 0;
}

size_t ReadFileVtk(const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, const std::string name_data,
	vtkDataArray*& data, const bool is_print/*=false*/) {

	/*„тение исходного файла и запись в vtkUnstructuredGrid*/

	vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
		vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader_vtk->ReadAllScalarsOn();
	reader_vtk->SetFileName(name_file_vtk.c_str());
	reader_vtk->Update();

	if (reader_vtk->IsFileUnstructuredGrid()) {
		unstructuredgrid = reader_vtk->GetUnstructuredGridOutput();
		unstructuredgrid->Modified();
	}
	else {
		std::cout << "Error read file\n";
		std::cout << "file_vtk is not UnstructuredGrid\n";
		return 1;
	}

	data = unstructuredgrid->GetCellData()->GetScalars(name_data.c_str());

	if (is_print) {
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfPoints() << " points.\n";
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfCells() << " cells.\n";

	}

	reader_vtk->GetOutput()->GlobalReleaseDataFlagOn();
	return 0;
}

size_t WriteFileSolutionEnergy1D(const std::string name_file_out, std::vector<int>& id_cells,
	vtkDataArray*& energy, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {

	std::ofstream ofile;
	ofile.open(name_file_out);
	if (!ofile.is_open()) {
		std::cout << "Error open file\n";
		std::cout << "file_graph is not opened for writing\n";
		return 1;
	}

	
	Vector3 center;
	for (auto cur_id: id_cells){
		CenterOfTetra(cur_id, unstructured_grid, center);
		ofile << center.norm() << " " << energy->GetTuple1(cur_id) << '\n';
	}

	ofile.close();

	return 0;
}

int GetFunction(const std::string name_file_vtk, const std::string name_file_out) {
	
	vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	vtkDataArray* Energy;

	ReadFileVtk(name_file_vtk, unstructured_grid, "scalars", Energy, true);

	std::vector<int> id_cell;

	GetIdTrace(unstructured_grid, id_cell);

	WriteFileSolutionEnergy1D(name_file_out, id_cell, Energy, unstructured_grid);

	return 0;
}