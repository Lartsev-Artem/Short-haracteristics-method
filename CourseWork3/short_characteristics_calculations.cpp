#include"short_characteristics_calculations.h"

Type BoundaryFunction(const Vector3 x) {
	return 0;
}

Type CalculateIllumeOnInnerFace(const int num_in_face, const int global_num_in_face, const Eigen::Matrix4d& vertex_tetra,
	const Vector3& x, const Vector3& x0, const std::map<int, Vector3>& nodes_value,
	const Matrix3& straight_face, const Matrix3& inclined_face, const Matrix3& transform_matrix, const Vector3& start_point_plane_coord) {
	Type I_x0;
	if (global_num_in_face == -1) {
		/*Граничные условия*/
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
		else {//наклонная плоскость

			FromGlobalToLocalTetra(vertex_tetra, x0, x0_local);
			Vector3 local_plane_x0;
			FromTetraToPlane(transform_matrix, start_point_plane_coord, x0_local, local_plane_x0);

			Vector3 coef = GetInterpolationCoef(inclined_face, nodes_value.find(global_num_in_face)->second);
			I_x0 = local_plane_x0[0] * coef[0] + local_plane_x0[1] * coef[1] + coef[2];
		}
	}
	return I_x0;
}

int CalculateNodeValue(const int num_cur_cell, vtkCell* cur_cell, const int num_cur_out_face, const int* face_state,
	const Vector3& direction,
	std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate,
	const int num_node, const Eigen::Matrix4d& vertex_tetra, Vector3& x,
	const Matrix3& straight_face, const Matrix3& inclined_face, const Matrix3& transform_matrix, const Vector3& start_point_plane_coord) {

	Vector3 x0;

	for (size_t num_in_face = 0; num_in_face < 4; ++num_in_face) {
		if (!face_state[num_in_face]) continue;  // обрабатываем только входные грани

		//if (direction.dot(normals.row(num_in_face)) < eps) {
		//	std::cout << "direction.normal==0, SOS!!!\n"; continue; // плоскость параллельна направлению
		//}

		IntersectionWithPlane(cur_cell->GetFace(num_in_face), x, direction, x0);

		if (InTriangle(cur_cell->GetFace(num_in_face), x0)) {

			//cout << all_pairs_face[num_cur_cell * 4 + num_in_face] << "  " << num_cur_cell * 4 + num_in_face << '\n';
			int global_num_in_face = Min(all_pairs_face[num_cur_cell * 4 + num_in_face], num_cur_cell * 4 + num_in_face);
			Type s = (x - x0).norm();

			// значение на входящей грани
			Type I_x0 = CalculateIllumeOnInnerFace(num_in_face, global_num_in_face, vertex_tetra, x, x0, nodes_value,
				straight_face, inclined_face, transform_matrix, start_point_plane_coord);

			int global_num_out_face = Min(all_pairs_face[num_cur_cell * 4 + num_cur_out_face], num_cur_cell * 4 + num_cur_out_face);
			if (global_num_out_face == -1)global_num_out_face = num_cur_cell * 4 + num_cur_out_face;

			nodes_value.find(global_num_out_face)->second[num_node] =
				GetIllum(num_cur_cell, s, I_x0, density, absorp_coef, rad_en_loose_rate);
			break;
		}

	}//for num_in_face

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

int ChooseDirection(const size_t number_direction, const vtkSmartPointer<vtkUnstructuredGrid>& sphere_direction_grid, Type* main_direction) {
	/*choosing the direction along the vertices of the sphere*/
	sphere_direction_grid->GetPoint(number_direction, main_direction);

	return 0;
}

size_t IntersectionWithPlane(vtkCell* face, const Vector3& start_point, const Vector3& direction, Vector3& result) {

	//вершины треугольника

	Type A[3];
	Type B[3];
	Type C[3];
	face->GetPoints()->GetPoint(0, A);
	face->GetPoints()->GetPoint(1, B);
	face->GetPoints()->GetPoint(2, C);


	Type a, b, c, d;  // параметры уравнения плоскости
	Type t;

	a = A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]);
	b = A[0] * (C[2] - B[2]) + B[0] * (A[2] - C[2]) + C[0] * (B[2] - A[2]);
	c = A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]);
	d = A[0] * (C[1] * B[2] - B[1] * C[2]) + B[0] * (A[1] * C[2] - C[1] * A[2]) + C[0] * (B[1] * A[2] - A[1] * B[2]);

	t = -(a * start_point[0] + b * start_point[1] + c * start_point[2] + d) / (a * direction[0] + b * direction[1] + c * direction[2]);

	for (size_t i = 0; i < 3; ++i)
		result[i] = (direction[i] * t + start_point[i]);  // точка пересечения луча  (start->direction) с плоскостью!!! face

	return 0;
}

size_t SetBasis(const Type* start_point, Vector3& normal, Matrix3& basis) {
	/*по начальной точке и нормале строит локальный базис картинной плоскости (vec1, vec2).
	  нормаль дана. задаем один вектор произвольно(ортагонально normal). третий вектор из векторного произведения*/
	Vector3 vec_1;
	Vector3 vec_2;

	if (abs(normal[1]) < 1e-20) {
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

	basis.col(1) = vec_1;
	basis.col(2) = vec_2;
	basis.col(3) = normal;

	return 0;
}
size_t Make2dPoint(const Type* start, const Matrix3& local_basis, const Type* point, Vector3& new_point) {



	for (size_t i = 0; i < 3; i++)
		new_point[i] = 0;

	//перевод 3d точки в 2d (в локальном базисе {start, local_basis}) 
	for (size_t k = 0; k < 3; k++) {
		new_point[0] += (point[k] - start[k]) * local_basis(0,k);
		new_point[1] += (point[k] - start[k]) * local_basis(1,k);
	}
	return 0;
}
bool InTriangle(vtkCell* face, const Eigen::Vector3d& XX) {
	/*face --- треугольник, X --- точка для проверки*/

	// вершины треугольника
	Type AA[3];
	Type BB[3];
	Type CC[3];
	face->GetPoints()->GetPoint(0, AA);
	face->GetPoints()->GetPoint(1, BB);
	face->GetPoints()->GetPoint(2, CC);

	
	Vector3 A, B, C, X;
	{
		Type Xx[3] = { XX[1],XX[2],XX[3] };
		Vector3 n;
		Matrix3 basis;
		NormalToFace(face, n);
		SetBasis(AA, n, basis);
		Make2dPoint(AA, basis, AA, A);
		Make2dPoint(AA, basis, BB, B);
		Make2dPoint(AA, basis, CC, C);
		Make2dPoint(AA, basis, Xx, X);
	}

	// линейная алгебра
	Type r1 = (A[0] - X[0]) * (B[1] - A[1]) - (B[0] - A[0]) * (A[1] - X[1]);
	Type r2 = (B[0] - X[0]) * (C[1] - B[1]) - (C[0] - B[0]) * (B[1] - X[1]);
	Type r3 = (C[0] - X[0]) * (A[1] - C[1]) - (A[0] - C[0]) * (C[1] - X[1]);

	if (r1 < 0 && r2 < 0 && r3 < 0)
		return true;
	else if (r1 > 0 && r2 > 0 && r3 > 0)
		return true;
	else return false;
}
Type IntegarteDirection(const vector<Type>& squares, const Type scuare_surface) {
	Type res = 0;
	int n = squares.size();
	for (size_t i = 0; i < n; ++i){
		res += 1 * squares[i] ;
	}
	return res / scuare_surface;
}
Type IntegarteDirection(const int num_cell, const vector<Type>& Illum, const vector<Type>& squares, const Type scuare_surface) {
	Type res = 0;
	int n = squares.size();
	
	for (size_t i = 0; i < n; i++) {
		res += Illum[n * i + num_cell] * squares[i];
		}

	return res / scuare_surface;
}
int InitGlobalValue(Vector3& start_point_plane_coord, Matrix3& transform_matrix, Matrix3& inverse_transform_matrix,
	Matrix3&	straight_face, Matrix3& inclined_face) {
	// 3 узла интерполяции
		{
			straight_face << 1. / 6, 1. / 6, 1,
				2. / 3, 1. / 6, 1,
				1. / 6, 2. / 3, 1;
		}

		// 3 узла интерполяции на наклонной плоскости
		{
			inclined_face <<
				0, sqrt(2. / 3), 1,
				sqrt(2) / 4, 1. / (2 * sqrt(6)), 1,
				-sqrt(2) / 4, 1. / (2 * sqrt(6)), 1;
		}

		//Матрицы перехода из стандартного тетраэдра в координаты наклонной плоскости 
		{ transform_matrix <<
			-1. / sqrt(2), 1. / sqrt(2), 0,
			-1. / sqrt(6), -1. / sqrt(6), sqrt(2. / 3),
			1. / sqrt(3), 1. / sqrt(3), 1. / sqrt(3);
		}

		//Матрицы перехода из наклонной плоскости в  координаты стандартного тетраэдра
		{
			inverse_transform_matrix <<
				-1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
				1. / sqrt(2), -1. / sqrt(6), 1. / sqrt(3),
				0, sqrt(2. / 3), 1. / sqrt(3);
		}

		// Начало координата плоскости
		start_point_plane_coord << 0.5, 0.5, 0;
	return 0;
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

	return count_unique_face;
}

int FromSphericalToDecart(const int number_cur, const vector<Type>& all_directions, Vector3& direction) {
	Type theta = all_directions[number_cur];
	Type fi = all_directions[all_directions.size() / 2 + number_cur];

	direction[0] = sin(theta) * cos(fi);
	direction[1] = sin(theta) * sin(fi);
	direction[2] = cos(theta);
	return 0;
}
int FromGlobalToLocalTetra(const Eigen::Matrix4d& vertex_tetra, const Eigen::Vector3d& global_coord, Eigen::Vector3d& local_coord)
{
	// vertex_tetra -> [X;Y;Z;1]
	// возможно надо будет использовать transpose из-за инициализации матрицы перехода

	Eigen::Matrix4d vertex_tetra_inverse = vertex_tetra.inverse();
	// сразу с транспонированием
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
	// написать в ручную т.к. преобразования известны, и 4я строка постоянна и меняться не должна

	Type eta4 = 1 - local_coord[0] - local_coord[1] - local_coord[2];
	for (int i = 0; i < 3; ++i)
	{
		global_coord[i] = 0;
		for (int j = 0; j < 3; ++j)
			global_coord[i] += vertex_tetra(i, j) * local_coord[j];
		global_coord[i] += vertex_tetra(i, 3) * eta4;
	}

	//global_coord = vertex_tetra.inverse() * local_coord;

	return 0;
}

int FromTetraToPlane(const Eigen::Matrix3d& transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& tetra_coord, Eigen::Vector3d& plane_coord) {
	plane_coord = transform_matrix * (tetra_coord - start_point);
	return 0;
}
int FromPlaneToTetra(const Eigen::Matrix3d& inverse_transform_matrix, const Eigen::Vector3d& start_point, const Eigen::Vector3d& plane_coord, 
	Eigen::Vector3d& tetra_coord) {
	tetra_coord = inverse_transform_matrix * plane_coord + start_point;
	return 0;
}

Type GetIllum(const int cur_id, const Type s, const Type I_node_prev,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate) {

	Type Ie = 10; // rad_en_loose_rate->GetTuple1(cur_id);
	Type k = 1;


	return Ie + (I_node_prev - Ie) * exp(-s * k);
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

int GetInterpolationCoef(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value, Eigen::Vector3d& coef) {
	//interpolation_nodes --- постояные узлы интерполяции (формат (в координатах стандартного тетраэдра) {x_i, y_i, 1})

	coef = interpolation_nodes.partialPivLu().solve(function_value);

	return 0;
}
Vector3 GetInterpolationCoef(const Eigen::Matrix3d& interpolation_nodes, const Eigen::Vector3d& function_value) {
	//interpolation_nodes --- постояные узлы интерполяции (формат (в координатах стандартного тетраэдра) {x_i, y_i, 1})
	return interpolation_nodes.partialPivLu().solve(function_value);
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

Type GetValueInCenterCell(const int num_cell, vtkCell* cur_cell, const Vector3 center, const Vector3 direction,
	const Eigen::Matrix4d& vertex_tetra, const std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate,
	const Matrix3& straight_face, const Matrix3& inclined_face, const Matrix3& transform_matrix, const Vector3& start_point_plane_coord) {
	/*Все грани должно быть определены*/
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


int GetNodesValues(const int num_cur_cell, vtkCell* cur_cell, const int num_cur_out_face, const int* face_state,
	const Vector3& direction, const Eigen::Matrix4d& vertex_tetra,
	std::map<int, Vector3>& nodes_value, const std::vector<int>& all_pairs_face,
	vtkDataArray* density, vtkDataArray* absorp_coef, vtkDataArray* rad_en_loose_rate,
	const Matrix3& straight_face, const Matrix3& inclined_face, const Matrix3& transform_matrix, const Matrix3& inverse_transform_matrix,
	const Vector3& start_point_plane_coord) {

	Vector3 x;

	Vector3 node;

	switch (num_cur_out_face)
	{
	case 0:
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = 0;
			node[1] = straight_face.row(num_node)[0];
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани

			CalculateNodeValue(num_cur_cell, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				density, absorp_coef, rad_en_loose_rate, num_node, vertex_tetra, x,
				straight_face, inclined_face, transform_matrix, start_point_plane_coord);
		}
		break;
	case 1:
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			node[0] = straight_face.row(num_node)[0];
			node[1] = 0;
			node[2] = straight_face.row(num_node)[1];
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани

			CalculateNodeValue(num_cur_cell, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				density, absorp_coef, rad_en_loose_rate, num_node, vertex_tetra, x,
				straight_face, inclined_face, transform_matrix, start_point_plane_coord);
		}
		break;
	case 2:
		for (size_t num_node = 0; num_node < 3; ++num_node) {			
		    node[0] = straight_face.row(num_node)[0];
			node[1] = straight_face.row(num_node)[1];
			node[2] = 0;
			FromLocalToGlobalTetra(vertex_tetra, node, x);  // x->координата узла на выходящей грани

			CalculateNodeValue(num_cur_cell, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				density, absorp_coef, rad_en_loose_rate, num_node, vertex_tetra, x,
				straight_face, inclined_face, transform_matrix, start_point_plane_coord);
		}
		break;
	case 3:
		for (size_t num_node = 0; num_node < 3; ++num_node) {		
			node[0] = inclined_face.row(num_node)[0];
			node[1] = inclined_face.row(num_node)[1];
			node[2] = 0;

			FromPlaneToTetra(inverse_transform_matrix, start_point_plane_coord, node, node);
			FromLocalToGlobalTetra(vertex_tetra, node, x);
		

			CalculateNodeValue(num_cur_cell, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				density, absorp_coef, rad_en_loose_rate, num_node, vertex_tetra, x,
				straight_face, inclined_face, transform_matrix, start_point_plane_coord);
		}
		break;
	default:
		std::cout << "Number face is not {0,1,2,3}????\n";
		break;
	}
	return 0;
	if (num_cur_out_face != 3) {
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			
			//Важен номер прямой грани. + значение 3 координаты 
			FromLocalToGlobalTetra(vertex_tetra, straight_face.row(num_node), x);  // x->координата узла на выходящей грани

			CalculateNodeValue(num_cur_cell, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				density, absorp_coef, rad_en_loose_rate, num_node, vertex_tetra, x,
				straight_face, inclined_face, transform_matrix, start_point_plane_coord);

		}//for num_node			
	}
	else {
		for (size_t num_node = 0; num_node < 3; ++num_node) {
			FromLocalToGlobalTetra(vertex_tetra, inclined_face.row(num_node), x);
			Vector3 local_tetra_x = x;
			FromTetraToPlane(transform_matrix, start_point_plane_coord, local_tetra_x, x);

			CalculateNodeValue(num_cur_cell, cur_cell, num_cur_out_face, face_state, direction, nodes_value, all_pairs_face,
				density, absorp_coef, rad_en_loose_rate, num_node, vertex_tetra, x,
				straight_face, inclined_face, transform_matrix, start_point_plane_coord);
		}
	}

	return 0;
}


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

int NormalsToCell(vtkCell* cell_tetra, Eigen::Matrix<Type, 4, 3>& normals) {

	Vector3 n;
	Type P0[3], P1[3], P2[3];
	for (int i = 0; i < 4; i++) {
		cell_tetra->GetFace(i)->GetPoints()->GetPoint(0, P0);
		cell_tetra->GetFace(i)->GetPoints()->GetPoint(1, P1);
		cell_tetra->GetFace(i)->GetPoints()->GetPoint(2, P2);

		Type a[3], b[3];
		for (size_t i = 0; i < 3; ++i) {
			a[i] = P1[i] - P0[i];
			b[i] = P2[i] - P0[i];
		}
		n[0] = a[1] * b[2] - a[2] * b[1];
		n[1] = -a[0] * b[2] + a[2] * b[0];
		n[2] = a[0] * b[1] - a[1] * b[0];

		n.normalize();

		for (size_t j = 0; j < 3; j++)
			normals(i, j) = n[j];
	}
	return 0;
}

int MakeFileCellData(const std::string name_file, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, bool rewrite_if_exist) {
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

	if (!rewrite_if_exist) {
		ifstream ifile;
		ifile.open(name_file);
		if (ifile.is_open()) {
			// файл уже существует
			int n;
			ifile >> n;
			if (n == N) {
				std::cout << "File CellData exist\n";
				return 0;
			}
		}
	}
	ofstream ofile;

	ofile.open(name_file);
	if (!ofile.is_open()) {
		std::cout << "error file CellData\n";
		return 1;
	}

	ofile << N << '\n';

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

size_t MakeEnergy(const vector<Type>& Illum, const vector<Type>& squares,const Type scuare_surface, vector<Type>& energy) {

	const int n = energy.size();

	for (size_t i = 0; i < n; ++i) {
		energy[i] = IntegarteDirection(i, Illum, squares, scuare_surface);
	}

	return 0;
}

int FindInAndOutFaces(const Vector3  direction, vtkCell* cur_cell_tetra, int* face_state) {
	//face_state  -0=> выходящая грань,  1=> входящая  face_state.size=4!!!  

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

int InitNodesValue(const std::vector<int>& all_pairs_face, std::map<int, Vector3>& nodes_value) {

	const int n = all_pairs_face.size();
	for (size_t i = 0; i < n; ++i)
	{
		if (i < all_pairs_face[i])
			nodes_value.insert({ i, Vector3(-2,-2,-2) });
	}
	return 0;
}

int SetVertexMatrix(const size_t number_cell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, Eigen::Matrix4d& vertex_tetra) {

	// 4 вершины треугольника(по столбцам и единицы в нижний строке)

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

	// 4 вершины треугольника(по столбцам и единицы в нижний строке)

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

	// 4 вершины треугольника(по столбцам и единицы в нижний строке)
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

size_t SortCellsGrid(Vector3 main_direction, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vector<int>& sorted_id_cell,
	const vector<Vector3>& centers_tetra) {

	Vector3 start_point = { -main_direction[0] * 5,-main_direction[1] * 5 ,-main_direction[2] * 5 };

	auto cmp{ [&unstructuredgrid,&main_direction,&start_point,&centers_tetra](const int left, const int right) {

		Vector3 left_point = centers_tetra[left];//{centers_tetra[left][0],centers_tetra[left][1],centers_tetra[left][2]};
		Vector3 right_point = centers_tetra[right];// { centers_tetra[right][0],centers_tetra[right][1],centers_tetra[right][2] };

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
		return L1 > L2;
	} };

	for (size_t i = 0; i < unstructuredgrid->GetNumberOfCells(); i++)
	{
		sorted_id_cell[i] = i;
	}

	sort(sorted_id_cell.begin(), sorted_id_cell.end(), cmp);

	return 0;
}
int FindCentersSpheres(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, const Vector3& local_center, vector<Vector3>& centers_tetra) {
	
	const int n = unstructured_grid->GetNumberOfCells();
	centers_tetra.resize(n);
	Eigen::Matrix4d vertex_tetra;
	

	for (size_t i = 0; i < n; ++i){
		SetVertexMatrix(i, unstructured_grid, vertex_tetra);
		FromLocalToGlobalTetra(vertex_tetra, local_center, centers_tetra[i]);

		cout << vertex_tetra << "\n\n\n" << centers_tetra[i] << '\n';
	}

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
	int N;
	ifile >> N;
	getline(ifile, str);

	centers_tetra.resize(N);

	for (int i = 0; i < N; ++i) {
		getline(ifile, str); getline(ifile, str); getline(ifile, str); getline(ifile, str); getline(ifile, str);
		ifile >> centers_tetra[i][0] >> centers_tetra[i][1] >> centers_tetra[i][2];
		getline(ifile, str);
	}

	ifile.close();
	return 0;
}
int ReadNextCellData(ifstream& ifs, Eigen::Matrix<Type, 4, 3>& normals, Vector3& center) {
	std::string str;
	getline(ifs, str); // порядок граней(сейчас в статичном варианте)
	for (size_t i = 0; i < 4; ++i) {
		for (size_t j = 0; j < 3; j++) {
			ifs >> normals(i, j);
		}
	}
	ifs >> center[0] >> center[1] >> center[2];
	getline(ifs, str);

	return 0;
}
int ReadCenters(const std::string name_file, vector<Vector3>& centers_tetra) {
	ifstream ifile;
	ifile.open(name_file);
	if (!ifile.is_open()) {
		std::cout << "File to read centers tetra  wasn't open\n";
		return 1;
	}

	std::string str;
	int N;
	ifile >> N;
	getline(ifile, str);

	centers_tetra.resize(N);

	for (int i = 0; i < N; ++i) {
		ifile >> centers_tetra[i][0] >> centers_tetra[i][1] >> centers_tetra[i][2];
		getline(ifile, str);
	}

	ifile.close();
	return 0;
}

size_t ReadFileVtk(const size_t class_file_vtk, const std::string name_file_vtk, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	vtkDataArray*& density, vtkDataArray*& absorp_coef, vtkDataArray*& rad_en_loose_rate, const bool is_print/*=false*/) {

	/*Чтение исходного файла и запись в vtkUnstructuredGrid*/

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

	switch (class_file_vtk) {
	case 0:
		density = NULL;
		absorp_coef = NULL;
		rad_en_loose_rate = NULL;
	case 1:
		density = unstructuredgrid->GetCellData()->GetScalars("alpha");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("alpha");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("Q");
		break;
	case 2:
		density = unstructuredgrid->GetCellData()->GetScalars("density");
		absorp_coef = unstructuredgrid->GetCellData()->GetScalars("absorp_coef");
		rad_en_loose_rate = unstructuredgrid->GetCellData()->GetScalars("radEnLooseRate");
		break;
	}

	if (is_print) {
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfPoints() << " points.\n";
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfCells() << " cells.\n";
		if (class_file_vtk) {
			std::cout << "density_Size: " << density->GetSize() << '\n';
			std::cout << "absorp_coef_Size: " << absorp_coef->GetSize() << '\n';
			std::cout << "Q_Size: " << rad_en_loose_rate->GetSize() << '\n';
		}
	}

	reader_vtk->GetOutput()->GlobalReleaseDataFlagOn();
	return 0;
}
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
size_t FromDecartToSphere(const vtkSmartPointer<vtkUnstructuredGrid>& sphere_direction_grid, vector<Type>& directions_all) {

	Type P[3];
	Type fi;
	Type theta;
	//  0---центр сферы
	int N = sphere_direction_grid->GetNumberOfPoints();
	for (size_t i = 1; i < N; i++) {
		sphere_direction_grid->GetPoint(i, P);

		FromDecartToSphere(P, fi, theta);
		directions_all[i] = theta;
		directions_all[(N - 1) + i] = fi;
	}
	return 0;
}

size_t ReadSphereDirectionVtk(const size_t class_file_vtk, const std::string name_file_sphere_direction, vector<Type>& directions_all) {

	/*Чтение исходного файла (vtk or txt) и запись в массив*/

	if (class_file_vtk == 0) {
		vtkSmartPointer<vtkUnstructuredGrid> unstructuredgrid =
			vtkSmartPointer<vtkUnstructuredGrid>::New();

		vtkSmartPointer<vtkGenericDataObjectReader> reader_vtk =
			vtkSmartPointer<vtkGenericDataObjectReader>::New();
		reader_vtk->ReadAllScalarsOn();
		reader_vtk->SetFileName(name_file_sphere_direction.c_str());
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

		const int N = unstructuredgrid->GetNumberOfPoints();
		directions_all.resize( N * 2);
		FromDecartToSphere(unstructuredgrid, directions_all);

		//sort(directions_all.begin(), directions_all.begin() + N);
		//sort(directions_all.begin() + N, directions_all.end());
	}
	else if (class_file_vtk == 1) {
		ifstream ifile;

		ifile.open(name_file_sphere_direction);
		if (!ifile.is_open()) {
			std::cout << "Error read file sphere direction\n";
			return 1;
		}
		int N = 0;
		ifile >> N;
		directions_all.resize(2 * N);

		int i = 0;
		while (!ifile.eof())
			ifile >> directions_all[i++] >> directions_all[i++];
		ifile.close();
	}
	else if (class_file_vtk == 2) {
		ifstream ifile;

		ifile.open(name_file_sphere_direction);
		if (!ifile.is_open()) {
			std::cout << "Error read file sphere direction\n";
			return 1;
		}
		int N = 0;
		ifile >> N;
		directions_all.resize(2 * N);

		int i = 0;
		Type P[3];
		Type theta;
		Type fi;
		while (!ifile.eof()) {
			ifile >> P[0] >> P[1] >> P[2];
			FromDecartToSphere(P, fi, theta);
			directions_all[i] = theta;
			directions_all[N + i++] = fi;
		}
		ifile.close();

		//sort(directions_all.begin(), directions_all.begin() + N);
		//sort(directions_all.begin() + N, directions_all.end());

	}
	else
		std::cout << "Error class_file_vtk\n";

	return 0;
}
size_t ReadSphereDirectionDecartToSpherical(const std::string name_file_sphere_direction, vector<Type>& directions_all, vector<Type>& squares, Type& square_surface) {

	ifstream ifile;

	ifile.open(name_file_sphere_direction);
	if (!ifile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}
	int N = 0;
	ifile >> N;
	directions_all.resize(2 * N);
	squares.resize(N);

	int i = 0;
	Type P[3];
	Type theta;
	Type fi;
	for (int i = 0; i < N; i++) {
		ifile >> squares[i];
		ifile >> P[0] >> P[1] >> P[2];
		FromDecartToSphere(P, fi, theta);
		directions_all[i] = theta;
		directions_all[N + i] = fi;
	}
	ifile >> square_surface;
	ifile.close();

	/**/

	//sort(directions_all.begin(), directions_all.begin() + N);
	//sort(directions_all.begin() + N, directions_all.end());

	return 0;
}

int MakeFileDirections(const std::string name_file_sphere_direction, vector<Type>& directions_all) {
	ofstream ofile;
	ofile.open(name_file_sphere_direction);
	if (!ofile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}

	int N = directions_all.size();
	ofile << N << '\n';

	/*sort(directions_all[0], directions_all[N / 2]);
	sort(directions_all[N / 2 + 1], directions_all[N - 1]);*/

	/*for (size_t i = 0; i < N; i+=2){
		ofile << directions_all[2 * i] << directions_all[2 * i + 1];
	}*/
	for (size_t i = 0; i < N; ++i) {
		ofile << fixed << setprecision(16) << directions_all[i] << '\n';
	}

	ofile.close();
	return 0;
}

size_t TransformFileDecartToSphere(const std::string name_file_sphere_direction, const std::string name_new_file_sphere_direction) {
	ifstream ifile;

	ifile.open(name_file_sphere_direction);
	if (!ifile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}

	vector<Type> directions_all;

	Type x, y, z;
	while (!ifile.eof()) {
		ifile >> x >> y >> z;
		directions_all.push_back(atan(sqrt(x * x + y * y) / z));
		directions_all.push_back(atan(y / x));
	}
	ifile.close();

	ofstream ofile;
	ofile.open(name_new_file_sphere_direction);

	ofile << directions_all.size() / 2 << '\n';
	for (size_t i = 0; i < directions_all.size(); i += 2) {
		ofile << directions_all[i] << ' ' << directions_all[i + 1] << '\n';
	}
	ofile.close();

	return 0;
}
size_t TransformNetgenToVtk(const std::string name_file_netgen, const std::string name_new_file_vtk) {

	ofstream ofile;
	ofile.open(name_new_file_vtk);
	if (!ofile.is_open()) {
		std::cout << "Error open file vtk\n";
		return 1;
	}

	ifstream ifile;
	ifile.open(name_file_netgen);
	if (!ifile.is_open()) {
		std::cout << "Error read file netgen\n";
		return 1;
	}

	ofile << "# vtk DataFile Version 3.0\n";
	ofile << "Elemensts Volumes\n";
	ofile << "ASCII\n";
	ofile << "DATASET UNSTRUCTURED_GRID\n";
	ofile << "POINTS ";

	int N;
	ifile >> N;
	ofile << N << " double\n";

	Type x, y, z;
	for (size_t i = 0; i < N; i++) {
		ifile >> x >> y >> z;
		ofile << x << ' ' << y << ' ' << z << '\n';
	}

	ifile >> N;
	ofile << "CELLS " << N << " " << N * 5 << '\n';

	int a, b, c, d;
	for (size_t i = 0; i < N; i++) {
		ifile >> a >> a >> b >> c >> d;
		ofile << 4 << ' ' << a - 1 << ' ' << b - 1 << ' ' << c - 1 << ' ' << d - 1 << '\n';
	}
	ofile << "CELL_TYPES " << N << '\n';
	for (size_t i = 0; i < N; i++)
		ofile << 10 << '\n';

	ofile.close();
	ifile.close();
	return 0;
}
size_t TransformNetgenToVtkSurface(const std::string name_file_netgen, const std::string name_new_file_vtk) {

	ofstream ofile;
	ofile.open(name_new_file_vtk);
	if (!ofile.is_open()) {
		std::cout << "Error open file vtk\n";
		return 1;
	}

	ifstream ifile;
	ifile.open(name_file_netgen);
	if (!ifile.is_open()) {
		std::cout << "Error read file netgen\n";
		return 1;
	}

	ofile << "# vtk DataFile Version 3.0\n";
	ofile << "Elemensts Volumes\n";
	ofile << "ASCII\n";
	ofile << "DATASET UNSTRUCTURED_GRID\n";
	ofile << "POINTS ";

	int N;
	ifile >> N;
	ofile << N << " double\n";

	Type x, y, z;
	for (size_t i = 0; i < N; i++) {
		ifile >> x >> y >> z;
		ofile << x << ' ' << y << ' ' << z << '\n';
	}

	ifile >> N;
	ofile << "CELLS " << N << " " << N * 4 << '\n';

	int a, b, c, d;
	for (size_t i = 0; i < N; i++) {
		ifile >> a >> a >> b >> c;
		ofile << 3 << ' ' << a - 1 << ' ' << b - 1 << ' ' << c - 1  << '\n';
	}
	ofile << "CELL_TYPES " << N << '\n';
	for (size_t i = 0; i < N; i++)
		ofile << 5 << '\n';

	ofile.close();
	ifile.close();
	return 0;
}

size_t PointIntersect(Eigen::Vector3d start_point, Eigen::Vector3d direction, Eigen::Vector3d& point, Eigen::Vector3d& result) {

	Type t = (-start_point.dot(direction) + direction.dot(point)) / (direction.dot(direction));

	result = direction * t + start_point; // точка пересечения луча  (start->direction) с плоскостью!!! face

	return 0;
}


size_t WriteFileSolution(const std::string NameFileOut, const std::vector<Type>& vector_energy, const vtkSmartPointer<vtkUnstructuredGrid>& UGrid) {

	int n = UGrid->GetNumberOfCells();

	vtkSmartPointer<vtkDoubleArray> IllumArray =
		vtkSmartPointer<vtkDoubleArray>::New();

	for (size_t i = 0; i < n; i++)
		IllumArray->InsertNextTuple1(vector_energy[i]);

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

int FindCentersSpheres(const vtkSmartPointer<vtkUnstructuredGrid>& unstructured_grid, vector<Vector3>& centers_tetra) {

	const int n = unstructured_grid->GetNumberOfCells();
	centers_tetra.resize(n);


	Type a[16];

	Eigen::MatrixXd Dx;
	Eigen::MatrixXd Dy;
	Eigen::MatrixXd Dz;
	Type detA;
	Vector3 center;
	Eigen::Vector4d xyz;

	for (size_t i = 0; i < n; ++i)
	{

		unstructured_grid->GetCell(i)->GetPoints()->GetPoint(0, a);
		a[3] = 1;
		unstructured_grid->GetCell(i)->GetPoints()->GetPoint(1, a + 4);
		a[7] = 1;
		unstructured_grid->GetCell(i)->GetPoints()->GetPoint(2, a + 8);
		a[11] = 1;
		unstructured_grid->GetCell(i)->GetPoints()->GetPoint(3, a + 12);
		a[15] = 1;

		Eigen::Matrix<Type, 4, 4> A(a);
		Dx = Dy = Dz = A;
		xyz[0] = A.col(0).dot(A.col(0)) - 1;
		xyz[1] = A.col(1).dot(A.col(1)) - 1;
		xyz[2] = A.col(2).dot(A.col(2)) - 1;
		xyz[3] = A.col(3).dot(A.col(3)) - 1;

		Dx.row(0) = xyz;
		Dy.row(0) = xyz;
		Dz.row(0) = xyz;

		Dy.row(1) = A.row(0);
		Dz.row(1) = A.row(0);
		Dz.row(2) = A.row(1);

		detA = 2 * A.determinant();

		center[0] = Dx.determinant() / detA;
		center[1] = -Dy.determinant() / detA;
		center[2] = Dz.determinant() / detA;
		centers_tetra[i] = center;

	}

	return 0;
}

int MakeFileDirectionsCenterTriangle(const std::string name_file_sphere_direction, const vtkSmartPointer<vtkUnstructuredGrid>& spehere_directions) {
	ofstream ofile;
	ofile.open(name_file_sphere_direction);
	if (!ofile.is_open()) {
		std::cout << "Error read file sphere direction\n";
		return 1;
	}

	int N = spehere_directions->GetNumberOfCells();
	ofile << N << '\n';
	
	auto MakeS{ [](Vector3& P0,Vector3& P1,Vector3& P2) {
		Type Sum = 0;
		Type a[3], b[3];
		for (size_t i = 0; i < 3; i++) {
			a[i] = P1[i] - P0[i];
			b[i] = P2[i] - P0[i];
		}

		Sum = pow(a[1] * b[2] - a[2] * b[1], 2) + pow(a[0] * b[2] - a[2] * b[0], 2) + pow(a[0] * b[1] - a[1] * b[0], 2);
		return 0.5 * sqrt(Sum);
} };


	Vector3 center;
	auto MakeCenter{ [&center](Vector3& P0,Vector3& P1,Vector3& P2) {
		
		Type a = (P0 - P1).norm();
		Type b = (P0 - P2).norm();
		Type c = (P2 - P1).norm();

		Type abc = a + b + c;
		center[0] = (a * P2[0] + b * P1[0] + c * P0[0]) / abc;
		center[1] = (a * P2[1] + b * P1[1] + c * P0[1]) / abc;
		center[2] = (a * P2[2] + b * P1[2] + c * P0[2]) / abc;

		return ;
} };


	//Type P0[3], P1[3], P2[3], P3[3];
	
	
	Type sum_squre = 0;
	Type cur_squre;
	for (size_t i = 0; i < N; ++i) {
		vtkSmartPointer<vtkIdList> idp = spehere_directions->GetCell(i)->GetPointIds();

		Vector3 P0(spehere_directions->GetPoint(idp->GetId(0)));
		Vector3 P1(spehere_directions->GetPoint(idp->GetId(1)));
		Vector3 P2(spehere_directions->GetPoint(idp->GetId(2)));

		cout << P0.norm() << '\n';


		/*spehere_directions->GetPoint(idp->GetId(1), P1);
		spehere_directions->GetPoint(idp->GetId(2), P2);*/
		cur_squre = MakeS(P0, P1, P2);
		sum_squre += cur_squre;
		MakeCenter(P0, P1, P2);

		ofile << fixed << setprecision(16) << cur_squre << '\n';
		ofile << fixed << setprecision(16) << center[0] << ' ' << center[1] << ' ' << center[2] << '\n';

		cout << center.norm() << '\n';
	}
	ofile << fixed << setprecision(16) << sum_squre;
	ofile.close();
	return 0;
}