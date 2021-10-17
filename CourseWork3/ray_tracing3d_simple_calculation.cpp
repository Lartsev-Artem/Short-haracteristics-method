#include "Header.h"
size_t CenterOfTetra(const int NumberCell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Type* PointInTetra);
size_t PointIntersect(Type* start_point, Type* direction, Type* point, Type* result);


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

size_t FromDecartToSphere(const vtkSmartPointer<vtkUnstructuredGrid>& sphere_direction_grid, vector<Type>& directions_all) {

	Type P[3];
	for (size_t i = 0; i < sphere_direction_grid->GetNumberOfPoints(); i++){
		sphere_direction_grid->GetPoint(i, P);
		directions_all[i] = atan(sqrt(P[0] * P[0] + P[1] * P[1]) / P[2]);
		directions_all[i + 1] = atan(P[1] / P[0]);
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

		directions_all.resize(unstructuredgrid->GetNumberOfPoints() * 2);
		FromDecartToSphere(unstructuredgrid, directions_all);
	}
	else {
		ifstream ifile;

		ifile.open(name_file_sphere_direction);
		if (!ifile.is_open()) {
			std::cout << "Error read file sphere direction\n";
			return 1;
		}
		int N = 0;
		ifile >> N;
		directions_all.resize(2*N);

		int i = 0;
		while (!ifile.eof())
			ifile >> directions_all[i++] >> directions_all[i++];
		ifile.close();	
	}

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
	for (size_t i = 0; i < directions_all.size(); i+=2){
		ofile << directions_all[i]<<' ' << directions_all[i + 1] << '\n';
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
	for (size_t i = 0; i < N; i++){
		ifile >> x >> y >> z;
		ofile << x << ' ' << y << ' ' << z << '\n';
	}

	ifile >> N;
	ofile << "CELLS " << N << " " << N * 5 << '\n';

	int a, b, c, d;
	for (size_t i = 0; i < N; i++) {
		ifile >> a >> a >> b >> c >> d;
		ofile << 4 << ' ' << a-1 << ' ' << b-1 << ' ' << c-1 << ' ' << d-1 << '\n';
	}
	ofile << "CELL_TYPES " << N << '\n';
	for (size_t i = 0; i < N; i++)
		ofile << 10 << '\n';

	ofile.close();
	ifile.close();
	return 0;
}

size_t SortCellsGrid(Type* main_direction, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, vector<int>& sorted_id_cell, vector<Eigen::Vector3d>& centers_tetra) {

	Type start_point[3] = { -main_direction[0] * 5,-main_direction[1] * 5 ,-main_direction[2] * 5 };

	auto cmp{ [&unstructuredgrid,&main_direction,&start_point,&centers_tetra](const int left, const int right) {

		Type left_point[3] = {centers_tetra[left][0],centers_tetra[left][1],centers_tetra[left][2]};
		Type right_point[3] = { centers_tetra[right][0],centers_tetra[right][1],centers_tetra[right][2] };
	
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

size_t PointIntersect(Type* start_point, Type* direction, Type* point,Type* result){
	
	Type a, b, c, d;  // параметры уравнения плоскости
	Type t;

	a = direction[0];
	b = direction[1];
	c = direction[2];
	d = -a * point[0] - b * point[1] - c * point[2];

	t = -(a * start_point[0] + b * start_point[1] + c * start_point[2] + d) / (a * direction[0] + b * direction[1] + c * direction[2]);

	for (size_t i = 0; i < 3; ++i)
		result[i] = (direction[i] * t + start_point[i]);  // точка пересечения луча  (start->direction) с плоскостью!!! face

	return 0;
}
size_t PointIntersect(Eigen::Vector3d start_point, Eigen::Vector3d direction, Eigen::Vector3d& point, Eigen::Vector3d& result) {

	Type t = (-start_point.dot(direction) + direction.dot(point)) / (direction.dot(direction));
	
	result = direction * t + start_point; // точка пересечения луча  (start->direction) с плоскостью!!! face

	return 0;
}

///*отладочная функция. Проверяет геометрическую длину пути луча внутри расчетной области*/
Type CheckLength(const std::vector<std::vector<Type>>& point) {
	if (point.size() == 0)
		return 0;
	Type sum = 0;
	Type lenght = 0;
	std::vector<Type> prev_point = point.front();
	
	for (auto el : point) {
		sum = 0;
		for (size_t i = 0; i < 3; ++i) {
			sum += pow(el[i] - prev_point[i], 2);
		}
		lenght += sqrt(sum);
		prev_point = el;
	//	std::cout << el[0] << ' ' << el[1] << ' ' << el[2] << '\n';
	}

	return lenght;
}

Type Distance(const Type* point1, const Type* point2) {
	Type s = 0;
	for (size_t i = 0; i < 3; ++i)
		s += pow((point2[i] - point1[i]), 2);
	return sqrt(s);
}
Type Distance(const Type point1_x, const Type point1_y, const Type point1_z,
	const Type point2_x, const Type point2_y, const Type point2_z) {

	Type s = pow(point1_x - point2_x, 2) + pow(point1_y - point2_y, 2) + pow(point1_z - point2_z, 2);
	return sqrt(s);
}

size_t IntersectionWithPlane(const std::vector<Type*>& face, const Type* start_point, const Type* direction, std::vector<Type>& result) {

	//вершины треугольника
	Type* A = face[0];
	Type* B = face[1];
	Type* C = face[2];

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

bool InTriangle(const std::vector<Type*>& face, const Type* X) {
	/*face --- треугольник, X --- точка для проверки*/

	// вершины треугольника
	Type* A = face[0];
	Type* B = face[1];
	Type* C = face[2];

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



size_t Make2dPoint(const Type* start, Type**& local_basis, const Type* point, Type* new_point) {

	for (size_t i = 0; i < 3; i++)
		new_point[i] = 0;

	//перевод 3d точки в 2d (в локальном базисе {start, local_basis}) 
	for (size_t k = 0; k < 3; k++) {
		new_point[0] += (point[k] - start[k]) * local_basis[0][k];
		new_point[1] += (point[k] - start[k]) * local_basis[1][k];
	}
	return 0;
}




Type MakeLength(Type* point1, Type* point2) {
	
	Type s = 0;  // безразмерная длина

	for (size_t i = 0; i < 3; ++i)
		s += (point2[i] - point1[i]) * (point2[i] - point1[i]);

	return sqrt(s);  // *388189 * pow(10, 5);  // числа --- переход к размерным параметрам
}

Type Norm(const Type* vec) {
	Type sum = 0;
	for (size_t i = 0; i < 3; ++i)
		sum += pow(vec[i], 2);
	return sqrt(sum);
}
size_t Normalize(Type* vec) {

	Type norm = Norm(vec);
	for (int i = 0; i < 3; i++)
		vec[i] /= norm;
	return 0;
}

Type Rosh(const Type* S, Type*& a, Type t) {

	// wolfram расчет
	Type max_potential = -1.82547;  // потенциал роша
	Type m = 0.8795180722891566;

	// математические расчёты
	Type x1 = S[0] + a[0] * t;
	Type x2 = pow(S[1] + a[1] * t, 2);
	Type x3 = pow(S[2] + a[2] * t, 2);

	return  -max_potential - (pow(x1 - m, 2) + x2 + x3) / 2 -
		m / sqrt(pow(x1 - 1, 2) + x2 + x3) -
		(1 - m) / sqrt(pow(x1, 2) + x2 + x3);
}

size_t SetBasis(const Type* start_point, const Type* normal, Type* vec_1, Type* vec_2) {
	/*по начальной точке и нормале строит локальный базис картинной плоскости (vec1, vec2).
	  нормаль дана. задаем один вектор произвольно(ортагонально normal). третий вектор из векторного произведения*/

	const Type end_point[3] = { 0.879518, 0, start_point[2] };  // параметры подобраны из геометрии задачи!!

	Type N[3];
	SetDirect(start_point, end_point, N);

	if (abs(N[1]) < pow(10, -20)) {
		vec_1[0] = 0;
		vec_1[1] = 1;
		vec_1[2] = 0;
	}
	else {
		vec_1[0] = 1;
		vec_1[2] = 0;
		vec_1[1] = -(N[0] * vec_1[0] + N[2] * vec_1[2]) / N[1];  //св-во скалярного произведения (N, vec1)==0
	}

	// правельная ориентация базиса плоскости
	if (normal[1] < 0)
		for (int i = 0; i < 3; ++i)
			vec_1[i] *= -1;

	// обычное векторное умножение. Eigen временно и не нужен!!!
	Eigen::Vector3d a(normal[0], normal[1], normal[2]);
	Eigen::Vector3d b(vec_1[0], vec_1[1], vec_1[2]);
	Eigen::Vector3d c = a.cross(b);

	for (size_t i = 0; i < 3; ++i)
		vec_2[i] = -c(i);

	Normalize(vec_1);
	Normalize(vec_2);

	return 0;
}

size_t SetDirect(const Type* start, const Type* end, Type* direct) {

	for (size_t i = 0; i < 3; ++i)
		direct[i] = (end[i] - start[i]);
	Normalize(direct);
	return 0;
}

size_t GetCenterTetra(const vtkCell* cell, Type* center) {

	Type A[3], B[3], C[3], D[3];
	vtkSmartPointer<vtkPoints> points = cell->Points;
	points->GetPoint(0, A);
	points->GetPoint(1, B);
	points->GetPoint(2, C);
	points->GetPoint(3, D);



	return 0;
}
size_t CenterOfTetra(const vtkCell* cell, Type* PointInTetra) {

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
	vtkSmartPointer<vtkPoints> points = cell->Points;
	points->GetPoint(0, P0);
	points->GetPoint(1, P1);
	points->GetPoint(2, P2);
	points->GetPoint(3, P3);

	Type Squr[4] = { MakeS(P1,P2,P3),MakeS(P0,P2,P3), MakeS(P0,P1,P3),MakeS(P0,P1,P2) };


	Type Sum = Squr[0] + Squr[1] + Squr[2] + Squr[3];
	for (size_t i = 0; i < 3; i++) {
		PointInTetra[i] = (Squr[0] * P0[i] + Squr[1] * P1[i] + Squr[2] * P2[i] + Squr[3] * P3[i]) / Sum;
	}
	return 0;
}

size_t CenterOfTetra(const int NumberCell, const vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Type* PointInTetra) {

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