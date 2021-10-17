#include"short_characteristics_calculations.h"

size_t SortCellsGrid(Vector3 main_direction, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vector<int>& sorted_id_cell,
	vector<Vector3>& centers_tetra) {

	Vector3 start_point = { -main_direction[0] * 5,-main_direction[1] * 5 ,-main_direction[2] * 5 };

	auto cmp{ [&unstructuredgrid,&main_direction,&start_point,&centers_tetra](const int left, const int right) {

		Vector3 left_point = {centers_tetra[left][0],centers_tetra[left][1],centers_tetra[left][2]};
		Vector3 right_point = { centers_tetra[right][0],centers_tetra[right][1],centers_tetra[right][2] };

		/*CenterOfTetra(left, unstructuredgrid, left_point);
		CenterOfTetra(right, unstructuredgrid, right_point);*/

		PointIntersect(start_point, main_direction, left_point, left_point);
		PointIntersect(start_point, main_direction, right_point, right_point);

		Type L1 = 0;
		Type L2 = 0;
		for (size_t i = 0; i < 3; ++i) {
			L1 += pow(start_point[i] - left_point[i], 2);
			L2 += pow(start_point[i] - right_point[i], 2);
		}
		return L1 < L2;
	} };

	for (size_t i = 0; i < unstructuredgrid->GetNumberOfCells(); i++)
	{
		sorted_id_cell[i] = i;
	}

	sort(sorted_id_cell.begin(), sorted_id_cell.end(), cmp);

	return 0;
}
int ReadCentersFromGeneralFile(const std::string name_file, vector<Vector3>& centers_tetra) {

	ifstream ifile;
	ifile.open(name_file);
	if (!ifile.is_open()) {
		std::cout << "File to read centers tetra  wasn't open\n";
		return 1;
	}

	std::string str;
	for (int i = 0; i < centers_tetra.size(); ++i) {
		getline(ifile, str); getline(ifile, str); getline(ifile, str); getline(ifile, str); getline(ifile, str);
		ifile >> centers_tetra[i][0] >> centers_tetra[i][1] >> centers_tetra[i][2];
		getline(ifile, str);
	}

	ifile.close();
	return 0;
}
Type CalculateIllumeOnInnerFace(const int num_in_face, const int global_num_in_face, const Eigen::Matrix4d& vertex_tetra,
	const Vector3& x, const Vector3& x0, const std::map<int, Vector3>& nodes_value, 
	const Matrix3& straight_face, const Matrix3& inclined_face, const Matrix3& transform_matrix, const Vector3& start_point_plane_coord) {
	Type I_x0;
	if (global_num_in_face == -1) {
		/*��������� �������*/
		I_x0 = BoundaryFunction(x0);
	}
	else {

		Vector3 x0_local;
		if (num_in_face != 3) {
			FromGlobalToLocalTetra(vertex_tetra, x0, x0_local);
			Vector3 coef = GetInterpolationCoef(straight_face, nodes_value.find(global_num_in_face)->second);
			switch (num_in_face) {
			case 0:
				I_x0 = x0_local[1] * coef[0] + x0_local[2] * coef[1] + coef[2];
				break;
			case 1:
				I_x0 = x0_local[0] * coef[0] + x0_local[2] * coef[1] + coef[2];
				break;
			case 2:
				I_x0 = x0_local[0] * coef[0] + x0_local[1] * coef[1] + coef[2];
				break;
			}
		}
		else {//��������� ���������

			FromGlobalToLocalTetra(vertex_tetra, x0, x0_local);
			Vector3 local_plane_x0;
			FromTetraToPlane(transform_matrix, start_point_plane_coord, x0_local, local_plane_x0);

			Vector3 coef = GetInterpolationCoef(inclined_face, nodes_value.find(global_num_in_face)->second);
			I_x0 = local_plane_x0[0] * coef[0] + local_plane_x0[1] * coef[1] + coef[2];
		}
	}
	return I_x0;
}

Type GetValueInCenterCell(const int num_cell, vtkCell* cur_cell, const Vector3 center, const Vector3 direction,
	const Eigen::Matrix4d& vertex_tetra, const std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate,
	const Matrix3& straight_face, const Matrix3& inclined_face, const Matrix3& transform_matrix, const Vector3& start_point_plane_coord) {
	/*��� ����� ������ ���� ����������*/
	Type value = 0;
	Vector3 x0;

	for (size_t i = 0; i < 4; i++) {

		IntersectionWithPlane(cur_cell->GetFace(i), center, direction, x0);
		if (InTriangle(cur_cell->GetFace(i), x0)) {

			Type s = (center - x0).norm();
			int global_num_in_face = Min(all_pairs_face[num_cell * 4 + i], num_cell * 4 + i);

			Type I_x0 = CalculateIllumeOnInnerFace(i, global_num_in_face, vertex_tetra, center, x0, nodes_value,
				straight_face, inclined_face,  transform_matrix, start_point_plane_coord);

			value = GetIllum(num_cell, s, I_x0, density, absorp_coef, rad_en_loose_rate);
			break;
		}
	}

	return value;
}

int CalculateNodeValue(const int num_cur_cell, vtkCell* cur_cell, const int num_cur_out_face, const int* face_state,
	const Vector3& direction,
	std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate,
	const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	const Matrix3& straight_face, const Matrix3& inclined_face, const Matrix3& transform_matrix, const Vector3& start_point_plane_coord) {

	Vector3 x0;

	for (size_t num_in_face = 0; num_in_face < 4; ++num_in_face) {
		if (!face_state[num_in_face]) continue;  // ������������ ������ ������� �����

		//if (direction.dot(normals.row(num_in_face)) < eps) {
		//	std::cout << "direction.normal==0, SOS!!!\n"; continue; // ��������� ����������� �����������
		//}

		IntersectionWithPlane(cur_cell->GetFace(num_in_face), x, direction, x0);

		if (InTriangle(cur_cell->GetFace(num_in_face), x0)) {

			int global_num_in_face = Min(all_pairs_face[num_cur_cell * 4 + num_in_face], num_cur_cell * 4 + num_in_face);
			Type s = (x - x0).norm();

			// �������� �� �������� �����
			Type I_x0 = CalculateIllumeOnInnerFace(num_in_face, global_num_in_face, vertex_tetra, x, x0, nodes_value,
				straight_face, inclined_face, transform_matrix, start_point_plane_coord);

			int global_num_out_face = Min(all_pairs_face[num_cur_cell * 4 + num_cur_out_face], num_cur_cell * 4 + num_cur_out_face);

			nodes_value.find(global_num_out_face)->second[num_node] =
				GetIllum(num_cur_cell, s, I_x0, density, absorp_coef, rad_en_loose_rate);
			break;
		}

	}//for num_in_face

	return 0;
}
int GetNodesValues(const int num_cur_cell, vtkCell* cur_cell, const int num_cur_out_face, const int* face_state,
	const Vector3& direction, const Eigen::Matrix4d& vertex_tetra,
	std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate,
	const Matrix3& straight_face, const Matrix3& inclined_face, const Matrix3& transform_matrix, const Vector3& start_point_plane_coord) {


	Vector3 x;

	if (num_cur_out_face != 3) {
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			FromLocalToGlobalTetra(vertex_tetra, straight_face.row(num_node), x);

			CalculateNodeValue(num_cur_cell, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				density, absorp_coef, rad_en_loose_rate, num_node, vertex_tetra, x,
				straight_face, inclined_face, transform_matrix, start_point_plane_coord);

		}//for num_node			
	}
	else {
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			FromLocalToGlobalTetra(vertex_tetra, straight_face.row(num_node), x);
			Vector3 local_tetra_x = x;
			FromTetraToPlane(transform_matrix, start_point_plane_coord, local_tetra_x, x);

			CalculateNodeValue(num_cur_cell, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				density, absorp_coef, rad_en_loose_rate, num_node, vertex_tetra, x,
				straight_face, inclined_face, transform_matrix, start_point_plane_coord);
		}
	}

	return 0;
}

int GetPointIdOppositeFace(vtkCell* cell_tetra, const size_t num_face) {

	vtkIdList* id_point_face = cell_tetra->GetFace(num_face)->GetPointIds();
	int id_face[3] = { id_point_face->GetId(0),id_point_face->GetId(1),id_point_face->GetId(2) };

	id_point_face = cell_tetra->GetPointIds();
	int id_cell[4] = { id_point_face->GetId(0),id_point_face->GetId(1),id_point_face->GetId(2),id_point_face->GetId(3) };

	bool is_equal = false;
	for (size_t i = 0; i < 4; i++)
	{
		bool is_equal = false;
		for (size_t j = 0; j < 3; j++)
		{
			if (id_cell[i] == id_face[j]) {
				is_equal = true;
				break;
			}

		}
		if (!is_equal)
			return id_cell[i];
	}

	cout << "Opposite point wasn't found??\n";
	return -1;
};
int NormalToFace(vtkCell* cell_face, Vector3& n) {

	vtkSmartPointer<vtkIdList> idp = cell_face->GetPointIds();

	Type P0[3], P1[3], P2[3];

	cell_face->GetPoints()->GetPoint(0, P0);
	cell_face->GetPoints()->GetPoint(1, P1);
	cell_face->GetPoints()->GetPoint(2, P2);

	Type a[3], b[3];
	for (size_t i = 0; i < 3; ++i) {
		a[i] = P1[i] - P0[i];
		b[i] = P2[i] - P0[i];
	}
	n[0] = a[1] * b[2] - a[2] * b[1];
	n[1] = -a[0] * b[2] + a[2] * b[0];
	n[2] = a[0] * b[1] - a[1] * b[0];

	n.normalize();
	return 0;
}
size_t CenterOfTetra(const int NumberCell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Vector3& PointInTetra) {

	auto MakeS{ [](Type* P0,Type* P1,Type* P2) {
		Type Sum = 0;
		Type a[3], b[3];
		for (size_t i = 0; i < 3; i++) {
			a[i] = P1[i] - P0[i];
			b[i] = P2[i] - P0[i];
		}

		Sum = pow(a[1] * b[2] - a[2] * b[1], 2) + pow(a[0] * b[2] - a[2] * b[0], 2) + pow(a[0] * b[1] - a[1] * b[0], 2);
		return 0.5 * sqrt(Sum);
} };

	Type P0[3], P1[3], P2[3], P3[3];
	vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetPointIds();

	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);
	unstructuredgrid->GetPoint(idp->GetId(2), P2);
	unstructuredgrid->GetPoint(idp->GetId(3), P3);

	Type Squr[4] = { MakeS(P1,P2,P3),MakeS(P0,P2,P3), MakeS(P0,P1,P3),MakeS(P0,P1,P2) };


	Type Sum = Squr[0] + Squr[1] + Squr[2] + Squr[3];
	for (size_t i = 0; i < 3; i++) {
		PointInTetra[i] = (Squr[0] * P0[i] + Squr[1] * P1[i] + Squr[2] * P2[i] + Squr[3] * P3[i]) / Sum;
	}
	return 0;
}
int MakeFileCellData(const std::string name_file, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid) {
	/*
	* FILE STRUCTURED
	* a b c d  --order vertex
	* n1  -- normal
	* n2
	* n3
	* n4
	* xx  --center tetra
	*/

	const int N = unstructured_grid->GetNumberOfCells();

	ofstream ofile;

	ofile.open(name_file);
	if (!ofile.is_open()) {
		std::cout << "error file CellData\n";
		return 1;
	}

	Vector3 vec3;

	for (size_t cur_cell = 0; cur_cell < N; ++cur_cell) {

		int i;
		for (size_t cur_face = 0; cur_face < 4; ++cur_face) {
			i = GetPointIdOppositeFace(unstructured_grid->GetCell(cur_cell), cur_face);
			ofile << i << ' ';
		}
		ofile << '\n';

		for (size_t i = 0; i < 4; i++) {
			NormalToFace(unstructured_grid->GetCell(cur_cell)->GetFace(i), vec3);
			ofile << fixed << setprecision(16) << vec3[0] << ' ' << vec3[1] << ' ' << vec3[2] << '\n';
		}

		CenterOfTetra(cur_cell, unstructured_grid, vec3);
		ofile << fixed << setprecision(16) << vec3[0] << ' ' << vec3[1] << ' ' << vec3[2] << '\n';
	}

	ofile.close();
	return 0;
}


int FindInAndOutFaces(const Vector3  direction, vtkCell* cur_cell_tetra, int* face_state) {
	//face_state  -0=> ��������� �����,  1=> ��������  face_state.size=4!!!  

	Vector3 normal;
	for (size_t i = 0; i < 4; ++i) {

		NormalToFace(cur_cell_tetra->GetFace(i), normal);

		if (normal.dot(direction) <= 0)
			face_state[i] = 1;
		else
			face_state[i] = 0;
	}

	return 0;
}

Type GetIllum(const int cur_id, const Type s, const Type I_node_prev,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate) {

	Type Ie = 10; // rad_en_loose_rate->GetTuple1(cur_id);
	Type k = 1;


	return Ie + (I_node_prev - Ie) * exp(-s * k);
}

int SetNodesValue(const std::vector<int>& all_pairs_face, std::map<int, Vector3>& nodes_value) {

	const int n = all_pairs_face.size();
	for (size_t i = 0; i < n; ++i)
	{
		if (i < all_pairs_face[i])
			nodes_value.insert({ i, Vector3() });
	}
	return 0;
}

int GetNumberNeighborFace(const int a, const int b, const int c, vtkCell* neighbor_cell) {

	vtkIdList* idc;

	int x, y, z;
	for (size_t i = 0; i < 4; i++)
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
}
int FindNeighborsPairFace(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, std::vector<int>& all_pairs_face) {

	int count_unique_face = 0;
	const int N = unstructured_grid->GetNumberOfCells();
	all_pairs_face.resize(N * 4);
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

			unstructured_grid->GetCellNeighbors(num_cell, idp, idc);

			if (idc->GetNumberOfIds() == 1) {
				int id_neighbor_cell = idc->GetId(0);
				int id_neighbor_face = GetNumberNeighborFace(id_a, id_b, id_c, unstructured_grid->GetCell(id_neighbor_cell));
				all_pairs_face[num_cell * 4 + num_face] = id_neighbor_cell * 4 + id_neighbor_face;
				all_pairs_face[id_neighbor_cell * 4 + id_neighbor_face] = num_cell * 4 + num_face;
			}
			else if (idc->GetNumberOfIds() == 0)
				all_pairs_face[num_cell * 4 + num_face] = -1; // ��������� ������
			else
				std::cout << "More than 1 neighbor????\n";
		}

	}

	return count_unique_face;
}

int FromSphericalToDecart(const int number_cur, const vector<Type>& all_directions, Vector3& direction) {
	Type theta = all_directions[number_cur];
	Type fi = all_directions[number_cur + 1];

	direction[0] = sin(theta) * cos(fi);
	direction[1] = sin(theta) * sin(fi);
	direction[2] = cos(theta);
	return 0;
}



int SetVertexMatrix(const size_t number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, Eigen::Matrix4d& vertex_tetra) {

	// 4 ������� ������������(�� �������� � ������� � ������ ������)

	for (size_t j = 0; j < 3; j++)
		vertex_tetra(j, 2) = unstructured_grid->GetCell(number_cell)->GetPoints()->GetPoint(2)[j];
	vertex_tetra(3, 2) = 1;

	for (size_t j = 0; j < 3; j++)
		vertex_tetra(j, 0) = unstructured_grid->GetCell(number_cell)->GetPoints()->GetPoint(0)[j];
	vertex_tetra(3, 0) = 1;

	for (size_t j = 0; j < 3; j++)
		vertex_tetra(j, 1) = unstructured_grid->GetCell(number_cell)->GetPoints()->GetPoint(1)[j];
	vertex_tetra(3, 1) = 1;

	for (size_t j = 0; j < 3; j++)
		vertex_tetra(j, 3) = unstructured_grid->GetCell(number_cell)->GetPoints()->GetPoint(3)[j];
	vertex_tetra(3, 3) = 1;

	return 0;
}
int SetVertexMatrix(vtkCell* cell, Eigen::Matrix4d& vertex_tetra) {

	// 4 ������� ������������(�� �������� � ������� � ������ ������)

	for (size_t j = 0; j < 3; j++)
		vertex_tetra(j, 2) = cell->GetPoints()->GetPoint(2)[j];
	vertex_tetra(3, 2) = 1;

	for (size_t j = 0; j < 3; j++)
		vertex_tetra(j, 0) = cell->GetPoints()->GetPoint(0)[j];
	vertex_tetra(3, 0) = 1;

	for (size_t j = 0; j < 3; j++)
		vertex_tetra(j, 1) = cell->GetPoints()->GetPoint(1)[j];
	vertex_tetra(3, 1) = 1;

	for (size_t j = 0; j < 3; j++)
		vertex_tetra(j, 3) = cell->GetPoints()->GetPoint(3)[j];
	vertex_tetra(3, 3) = 1;

	return 0;
}

int SetVertexMatrix(ifstream& ifs, const size_t number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, Eigen::Matrix4d& vertex_tetra) {

	// 4 ������� ������������(�� �������� � ������� � ������ ������)
	int num_vertex;
	for (size_t i = 0; i < 4; i++)
	{
		ifs >> num_vertex;
		for (size_t j = 0; j < 3; j++)
			vertex_tetra(j, i) = unstructured_grid->GetPoint(num_vertex)[j]; //unstructured_grid->GetCell(number_cell)->GetPoints()->GetPoint(num_vertex)[j];
		vertex_tetra(3, i) = 1;
	}
	return 0;
}

int ChooseDirection(const size_t number_direction, const vtkSmartPointer<vtkUnstructuredGrid>& sphere_direction_grid, Type* main_direction) {
	/*choosing the direction along the vertices of the sphere*/
	sphere_direction_grid->GetPoint(number_direction, main_direction);

	return 0;
}

int FromGlobalToLocalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& global_coord, Eigen::Vector3d& local_coord)
{
	// vertex_tetra -> [X;Y;Z;1]
	// �������� ���� ����� ������������ transpose ��-�� ������������� ������� ��������

	Eigen::Matrix4d vertex_tetra_inverse = vertex_tetra.inverse();
	// ����� � �����������������
	for (int i = 0; i < 3; ++i)
	{
		local_coord[i] = 0;
		for (int j = 0; j < 3; ++j)
			local_coord[i] += vertex_tetra_inverse(i, j) * global_coord[j];
		local_coord[i] += vertex_tetra_inverse(i, 3);
	}

	//local_coord = vertex_tetra * global_coord;
	return 0;
}
int FromLocalToGlobalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& local_coord, Eigen::Vector3d& global_coord) {
	// vertex_tetra -> [X,Y,Z,1]
	// �������� � ������ �.�. �������������� ��������, � 4� ������ ��������� � �������� �� ������

	Type eta4 = 1 - local_coord[0] - local_coord[1] - local_coord[2];
	for (int i = 0; i < 3; ++i)
	{
		global_coord[i] = 0;
		for (int j = 0; j < 3; ++j)
			global_coord[i] = vertex_tetra(i, j) * local_coord[j];
		global_coord[i] = vertex_tetra(i, 3) * eta4;
	}

	//global_coord = vertex_tetra.inverse() * local_coord;

	return 0;
}

int FromTetraToPlane(const Eigen::Matrix3d& transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& tetra_coord, Eigen::Vector3d& plane_coord) {
	plane_coord = transform_matrix * (tetra_coord - start_point);
	return 0;
}
int FromPlaneToTetra(const Eigen::Matrix3d& inverse_transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& plane_coord, Eigen::Vector3d& tetra_coord) {
	tetra_coord = inverse_transform_matrix * plane_coord + start_point;
	return 0;
}

int GetInterpolationCoef(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value, Eigen::Vector3d& coef) {
	//interpolation_nodes --- ��������� ���� ������������ (������ (� ����������� ������������ ���������) {x_i, y_i, 1})

	coef = interpolation_nodes.partialPivLu().solve(function_value);

	return 0;
}
Vector3 GetInterpolationCoef(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value) {
	//interpolation_nodes --- ��������� ���� ������������ (������ (� ����������� ������������ ���������) {x_i, y_i, 1})
	return interpolation_nodes.partialPivLu().solve(function_value);
}

int ReadNextCellData(ifstream& ifs, Eigen::Matrix<Type, 4, 3>& normals, Vector3& center) {
	std::string str;
	getline(ifs, str);
	for (size_t i = 0; i < 4; ++i) {
		for (size_t j = 0; j < 3; j++) {
			ifs >> normals(i, j);
		}
	}
	ifs >> center[0] >> center[1] >> center[2];
	getline(ifs, str);

	return 0;
}
size_t IntersectionWithPlane(vtkCell* face, const Vector3& start_point, const Vector3& direction, Vector3& result) {

	//������� ������������

	Type A[3];
	Type B[3];
	Type C[3];
	face->GetPoints()->GetPoint(0, A);
	face->GetPoints()->GetPoint(1, B);
	face->GetPoints()->GetPoint(2, C);


	Type a, b, c, d;  // ��������� ��������� ���������
	Type t;

	a = A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]);
	b = A[0] * (C[2] - B[2]) + B[0] * (A[2] - C[2]) + C[0] * (B[2] - A[2]);
	c = A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]);
	d = A[0] * (C[1] * B[2] - B[1] * C[2]) + B[0] * (A[1] * C[2] - C[1] * A[2]) + C[0] * (B[1] * A[2] - A[1] * B[2]);

	t = -(a * start_point[0] + b * start_point[1] + c * start_point[2] + d) / (a * direction[0] + b * direction[1] + c * direction[2]);

	for (size_t i = 0; i < 3; ++i)
		result[i] = (direction[i] * t + start_point[i]);  // ����� ����������� ����  (start->direction) � ����������!!! face

	return 0;
}
bool InTriangle(vtkCell* face, const Eigen::Vector3d& X) {
	/*face --- �����������, X --- ����� ��� ��������*/

	// ������� ������������
	Type A[3];
	Type B[3];
	Type C[3];
	face->GetPoints()->GetPoint(0, A);
	face->GetPoints()->GetPoint(1, B);
	face->GetPoints()->GetPoint(2, C);

	// �������� �������
	Type r1 = (A[0] - X[0]) * (B[1] - A[1]) - (B[0] - A[0]) * (A[1] - X[1]);
	Type r2 = (B[0] - X[0]) * (C[1] - B[1]) - (C[0] - B[0]) * (B[1] - X[1]);
	Type r3 = (C[0] - X[0]) * (A[1] - C[1]) - (A[0] - C[0]) * (C[1] - X[1]);

	if (r1 < 0 && r2 < 0 && r3 < 0)
		return true;
	else if (r1 > 0 && r2 > 0 && r3 > 0)
		return true;
	else return false;
}

Type BoundaryFunction(const Vector3 x) {
	return 0;
}

Type IntegarteDirection(const vector<Type>& directions) {
	Type res = 0;

	vector<Type> fi(directions.size() / 2);
	vector<Type> theta(fi.size());
	for (size_t i = 0; i < fi.size(); i++) {
		theta[i] = directions[i * 2];
		fi[i] = directions[i * 2 + 1];
	}

	sort(fi.begin(), fi.end());
	sort(theta.begin(), theta.end());

	for (size_t i = 0; i < fi.size() - 1; i++)
		for (size_t j = 0; j < theta.size() - 1; j++) {
			res += sin(theta[j]) * (theta[j + 1] - theta[j]) * (fi[i + 1] - fi[i]);
		}

	return res;
}