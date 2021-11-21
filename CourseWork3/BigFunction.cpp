#include "short_characteristics_main.h"
extern Vector3 start_point_plane_coord; // начало координат плоскости
extern Matrix3 transform_matrix;   // матрица перехода из базового тетраэдра в плоскость
extern Matrix3 inverse_transform_matrix;  // матрица перехода из плоскости в базовый тетраэдр

extern Matrix3	straight_face;  // 3 узла интерпол€ции
extern Matrix3 inclined_face;  // 3 узла интерпол€ции на наклонной плоскости

Type GetS(const int num_cell, const Vector3& direction, const std::vector<Type>& illum_old,
	const vector<Type>& directions, const vector<Type>& squares, const Type square_surface) {
	//num_cell equals x
	auto Gamma{ [](const Vector3& direction, const Vector3& direction2) {
	return (3. * (1 + pow(direction.dot(direction2),2))) / 4;
	} };


	Vector3 cur_direction;
	Type S = 0;
	const int N_dir = directions.size() / 2;
	const int N_cell = illum_old.size() / N_dir;

	for (int num_direction = 0; num_direction < N_dir; num_direction++)
	{
		FromSphericalToDecart(num_direction, directions, cur_direction);
		S += Gamma(cur_direction, direction) * illum_old[num_direction * N_cell + num_cell] * squares[num_direction];
	}
	return 4*PI*S / square_surface;
}

Type CurGetIllum(const int cur_id, const Vector3 x, const Type s, const Type I_node_prev, const Vector3& cur_direction,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate, const std::vector<Type>& illum_old,
	const vector<Type>& directions, const vector<Type>& squares, const Type square_surface) {

	Type Ie = 10 + GetS(cur_id, cur_direction, illum_old, directions, squares, square_surface);
	Type k = 10;
	if (x.norm() > 0.3) { Ie = 0; k = 1; }

	Type I;
	if (s > 1e-10)
		I = Ie * (1 - exp(-s * k)) + I_node_prev * exp(-s * k);
	else
		I = Ie * (1 + s * k) - I_node_prev * s * k;

	if (I < 0)
		I = 0;
	return I;
}

int CalculateNodeValue(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell,
	const int num_cur_out_face, const int* face_state,
	const Vector3& direction,
	std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face,
	const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate, const std::vector<Type>& illum_old,
	const vector<Type>& directions, const vector<Type>& squares, const Type square_surface) {

	Vector3 x0;

	for (size_t num_in_face = 0; num_in_face < 4; ++num_in_face) {
		if (!face_state[num_in_face]) continue;  // обрабатываем только входные грани


		IntersectionWithPlane(cur_cell->GetFace(num_in_face), x, direction, x0);

		if (InTriangle(num_cur_cell, unstructuredgrid, cur_cell, num_in_face, x0)) {

			int global_num_in_face = Min(all_pairs_face[num_cur_cell * 4 + num_in_face], num_cur_cell * 4 + num_in_face);

			Type s = (x - x0).norm();

			// значение на вход€щей грани
			Type I_x0 = CalculateIllumeOnInnerFace(num_in_face, global_num_in_face, vertex_tetra, x, x0, nodes_value,
				straight_face, inclined_face, transform_matrix, start_point_plane_coord);

			int id = num_cur_cell * 4 + num_cur_out_face;
			int global_num_out_face = Min(all_pairs_face[id], id);
			if (global_num_out_face == -1) global_num_out_face = id;


			Type I = CurGetIllum(num_cur_cell, x0, s, I_x0, direction, density, absorp_coef, rad_en_loose_rate, illum_old, directions, squares, square_surface);
			nodes_value.find(global_num_out_face)->second[num_node] = I;

			break;
		}

	}//for num_in_face

	return 0;
}

int GetNodes(const int num_cur_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vtkCell* cur_cell, const int num_cur_out_face,
	const Eigen::Matrix4d& vertex_tetra, const int* face_state, const Vector3& direction,
	std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate, const std::vector<Type>& illum_old,
	const vector<Type>& directions, const vector<Type>& squares, const Type square_surface) {

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

			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				num_node, vertex_tetra, x, density, absorp_coef, rad_en_loose_rate, illum_old, directions, squares, square_surface);
		}
		//cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 2://2->0
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = 0;
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выход€щей грани				
			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				num_node, vertex_tetra, x, density, absorp_coef, rad_en_loose_rate, illum_old, directions, squares, square_surface);
		}
		//cout << "Number: " << num_cur_out_face << "\n";
		break;
	case 0: //0->3
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = straight_face.row(num_node)[1];
			node[2] = 0;
			FromLocalToGlobalTetra(vertex_tetra, node, x);
			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				num_node, vertex_tetra, x, density, absorp_coef, rad_en_loose_rate, illum_old, directions, squares, square_surface);
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
			CalculateNodeValue(num_cur_cell, unstructuredgrid, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				num_node, vertex_tetra, x, density, absorp_coef, rad_en_loose_rate, illum_old, directions, squares, square_surface);
		}
		//	cout << "Number: " << num_cur_out_face << "\n";
		break;
	default:
		std::cout << "Number face is not {0,1,2,3}????\n";
		break;
	}

	return 0;
}

Type NormIllum(const std::vector<Type>& Illum, const std::vector<Type>& Illum2) {
	Type max = -1;
	Type buf;
	for (size_t i = 0; i < Illum.size(); i++)
	{
		buf = fabs(Illum[i] - Illum2[i]);
		if (buf > max)
			max = buf;
	}
	return max;
}

int Start(const std::string name_file_graph, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, 
	const std::vector<int>& all_pairs_face, std::map<int, Vector3>& nodes_value,
	const vector<Type>& directions, const vector<Type>& squares, const Type square_surface,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate,
	std::vector<Type>& Illum1) {


	std::vector<Type> Illum2(Illum1.size(), 0);
	vector<IntId> sorted_id_cell(unstructured_grid->GetNumberOfCells());
	int count = 0;

	do {
		Type _clock = -omp_get_wtime();
		{
			
			const int count_direction = directions.size() / 2;
			Eigen::Matrix4d vertex_tetra;
			Vector3 direction;

			/*---------------------------------- далее FOR по направлени€м----------------------------------*/
			for (int num_direction = 0; num_direction < count_direction; ++num_direction)
			{
				FromSphericalToDecart(num_direction, directions, direction);
				ReadGraph(name_file_graph + to_string(num_direction) + ".txt", sorted_id_cell);
				InitNodesValue(all_pairs_face, nodes_value);

				Vector3 x;
				Vector3 x0;
				int num_cell;

				/*---------------------------------- далее FOR по €чейкам----------------------------------*/
				for (int h = 0; h < unstructured_grid->GetNumberOfCells(); ++h) {
					num_cell = sorted_id_cell[h];

					SetVertexMatrix(num_cell, unstructured_grid, vertex_tetra);

					int face_state[4];  // 0=> выход€ща€ грань,  1=> вход€ща€   

					FindInAndOutFaces(direction, num_cell, unstructured_grid, face_state);

					for (size_t num_out_face = 0; num_out_face < 4; ++num_out_face) {
						if (!face_state[num_out_face])  // выход€щие грани
							GetNodes(num_cell, unstructured_grid, unstructured_grid->GetCell(num_cell), num_out_face, vertex_tetra, face_state, direction,
								nodes_value, all_pairs_face, density, absorp_coef, rad_en_loose_rate, Illum2, directions, squares, square_surface);
					}

					Vector3 center;
					CenterOfTetra(num_cell, unstructured_grid, center);
					Type I_k_dir = GetValueInCenterCell(num_cell, unstructured_grid, unstructured_grid->GetCell(num_cell), center, direction, vertex_tetra,
						nodes_value, all_pairs_face,
						density, absorp_coef, rad_en_loose_rate,
						straight_face, inclined_face, transform_matrix, start_point_plane_coord);

					Illum1[num_direction * unstructured_grid->GetNumberOfCells() + num_cell] = I_k_dir;
				}
				/*---------------------------------- конец FOR по €чейкам----------------------------------*/

				//std::cout << "End direction number: " << num_direction << '\n';
			}
			/*---------------------------------- конец FOR по направлени€м----------------------------------*/
		}

		std::swap(Illum1, Illum2);
		_clock += omp_get_wtime();
		count++;

		std::cout << "Time of iter: " << _clock << '\n';
		std::cout << "End iter_count number: " << count << '\n';
	} while (NormIllum(Illum1, Illum2) > 5e-2 && count < 100);
	
	std::cout << "Count:= " << count << '\n';
	std::cout << "Error:= " << NormIllum(Illum1, Illum2) << '\n';

	return 0;
}