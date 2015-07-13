#include <stdlib.h>
#include <stdio.h>
// #include <../../../prmt_sintactic_addition/prmt_sintactic_addition.h>

#include "../../../calculation_core/src/blocks/general/domain/domain.h"
#include "../../../calculation_core/src/blocks/general/laplacian/scalar/laplacian_scalar.h"
#include "../../../calculation_core/src/blocks/general/source/scalar/source_scalar.h"
#include "../../../calculation_core/src/blocks/general/boundary_value/boundary_value.h"
#include "../../../calculation_core/src/blocks/general/assembler/assembler.h"
#include "../../../calculation_core/src/blocks/general/system_linear_algebraic_equations/system_linear_algebraic_equations.h"
#include "../../../calculation_core/src/blocks/general/additional_tools/trivial_prepare_system_equations/trivial_prepare_system_equations.h"
#include "../../../calculation_core/src/blocks/general/additional_tools/apply_boundary_value/scalar/apply_boundary_value_scalar.h"
#include "../../../calculation_core/src/blocks/general/geometric_tools/geometric_tools.h"
#include "../../../calculation_core/src/blocks/special/heat_conduction_problem_tools/heat_conduction_problem_tools.h"
//#include "gradient.h"

// #include "../../../calculation_core/src/blocks/special/problem_on_cell/domain_looper/domain_looper.h"
// #include "../../../calculation_core/src/blocks/special/problem_on_cell/black_on_white_substituter/black_on_white_substituter.h"
#include "../../../calculation_core/src/blocks/special/problem_on_cell/source/scalar/source_scalar.h"
#include "../../../calculation_core/src/blocks/special/problem_on_cell/prepare_system_equations/prepare_system_equations.h"
#include "../../../calculation_core/src/blocks/special/problem_on_cell/prepare_system_equations_with_cubic_grid/prepare_system_equations_with_cubic_grid.h"
#include "../../../calculation_core/src/blocks/special/problem_on_cell/system_linear_algebraic_equations/system_linear_algebraic_equations.h"
#include "../../../calculation_core/src/blocks/special/problem_on_cell/calculate_meta_coefficients/calculate_meta_coefficients.h"
#include "../../../calculation_core/src/blocks/special/problem_on_cell/assembler/assembler.h"
// #include "../../../calculation_core/src/blocks/special/problem_on_cell/domain_looper_trivial/domain_looper_trivial.h"

#include "../../../calculation_core/src/blocks/general/laplacian/vector/laplacian_vector.h"
#include "../../../calculation_core/src/blocks/general/source/vector/source_vector.h"
#include "../../../calculation_core/src/blocks/general/additional_tools/apply_boundary_value/vector/apply_boundary_value_vector.h"
#include "../../../calculation_core/src/blocks/special/elastic_problem_tools/elastic_problem_tools.h"

#include "../../../calculation_core/src/blocks/special/problem_on_cell/source/vector/source_vector.h"
#include "nikola_elastic_problem_tools.h"


#include "../../../calculation_core/src/blocks/special/nikola_problem/source/scalar/source_scalar.h"
#include "../../../calculation_core/src/blocks/special/nikola_problem/source/vector/source_vector.h"

#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/grid/grid_reordering.h>

#include <deal.II/grid/grid_tools.h>						//для gmsh
#include <deal.II/grid/grid_in.h>						//для gmsh

#include "read_file_to_array.h"


namespace NEW_STRUCT_ELASTIC_T2_1
{
	void Add_to_vector(dealii::Vector<dbl> &Vec, const double number)						//Позволяет добавить к решению константу
	{
		for (int i=0; i<Vec.size(); ++i)	
		{
			Vec(i) += number;
		}
	}

	double Integral(const dealii::Vector<dbl> &Vec, Domain<2> &domain, dealii::FEValues<2> &fe_values, cst &n_q_points, cst &dofs_per_cell)						//Позволяет добавить к решению константу
	{
		vec<dealii::types::global_dof_index> local_dof_indices (dofs_per_cell);
		double integral_U_z = 0.0;
		double Square = 0.0;
//		std::cout << "\n";
//		std::cout << "n_q_points = " << n_q_points << "\n";											//4
//		std::cout << "dofs_per_cell = " << dofs_per_cell << "\n";									//8
//		std::cout << "local_dof_indices.size() = " << local_dof_indices.size() << "\n";				//8
//		std::cout << "(before cycle) Vec(10) = " << Vec(10) << "\n";
		for (auto cell = domain.dof_handler.begin_active(); cell != domain.dof_handler.end(); ++cell)
		{
			fe_values .reinit (cell);
			cell->get_dof_indices(local_dof_indices);
			for(int i = 0; i < n_q_points; ++i)
			{
				integral_U_z += Vec(local_dof_indices[2*i] / 2)*fe_values.JxW(i);
//				integral_U_z += Vec(local_dof_indices[i])*fe_values.JxW(i);
				Square += fe_values.JxW(i);
//				std::cout << "Vec(" << local_dof_indices[2*i] / 2 << ") = " << Vec(local_dof_indices[2*i] / 2) << "\t\t";
//				std::cout << "fe_values.JxW(i) = " << fe_values.JxW(i) << "\t\t";
//				std::cout << "Vec(" << local_dof_indices[i] << ")*fe_values.JxW(i) = " << Vec(local_dof_indices[i])*fe_values.JxW(i) << "\n";
			}
		}
//		std::cout << "0 / 2 = " << 0 / 2 << "  0\n";
//		std::cout << "1 / 2 = " << 1 / 2 << "  0\n";
//		std::cout << "2 / 2 = " << 2 / 2 << "  1\n";
//		std::cout << "3 / 2 = " << 3 / 2 << "  1\n";
//		std::cout << "4 / 2 = " << 4 / 2 << "  2\n";
//		std::cout << "5 / 2 = " << 5 / 2 << "  2\n";
//		std::cout << "6 / 2 = " << 6 / 2 << "  3\n";
//		std::cout << "7 / 2 = " << 7 / 2 << "  3\n";
//		std::cout << "8 / 2 = " << 8 / 2 << "  4\n";
//		std::cout << "9 / 2 = " << 9 / 2 << "  4\n";
/*
		std::cout << "Vec(10) = " << Vec(10) << "\n";
		std::cout << "Vec.size() = " << Vec.size() << "\n";
		std::cout << "Integral(U_z) = " << integral_U_z << "\n";
		std::cout << "Square = " << Square << "\n";
*/
//		cout << ".";
		return integral_U_z;
	}

	double IntegralOfCoordinateX(const dealii::Vector<dbl> &Vec, Domain<2> &domain, dealii::FEValues<2> &fe_values, cst &n_q_points, cst &dofs_per_cell, const double &X, const double &Yright, const double &Yleft, const int &N)						//Высчитывает интеграл по поверхности вдоль координаты
	{
		//N			-	количество отрезков на которые делится участок вычисления интеграла
		double Integral = 0.0;
		double DeltaY = (Yright - Yleft) / N;		//длинна отрезка
		double Fun = 0.0;							//Значение функции в середине отрезка
		dealii::Point<2, double> p(X, Yleft);

		for (int i = 0; i < N; ++i)
		{
			p(1) = Yleft + DeltaY / 2 + i*DeltaY;
			Fun = HCPTools::print_function_elastic<2>(Vec, domain, p);
			Integral += Fun * DeltaY;
//			std::cout << "Fun = " << Fun << "\n";
		}
//		std::cout << "Yright - Yleft = " << Yright - Yleft << "\n";		
//		std::cout << "DeltaY = " << DeltaY << "\n";
//		std::cout << "Integral = " << Integral << "\n";

		return Integral;
	}

	double IntegralOfCoordinateY(const dealii::Vector<dbl> &Vec, Domain<2> &domain, dealii::FEValues<2> &fe_values, cst &n_q_points, cst &dofs_per_cell, const double &Y, const double &Xup, const double &Xdown, const int &N)						//Высчитывает интеграл по поверхности вдоль координаты
	{
		//N			-	количество отрезков на которые делится участок вычисления интеграла
		double Integral = 0.0;
		double DeltaX = (Xup - Xdown) / N;		//длинна отрезка
		double Fun = 0.0;							//Значение функции в середине отрезка
		dealii::Point<2, double> p(Xdown, Y);

		for (int i = 0; i < N; ++i)
		{
			p(0) = Xdown + DeltaX / 2 + i*DeltaX;
			Fun = HCPTools::print_function_elastic<2>(Vec, domain, p);
			Integral += Fun * DeltaX;
//			std::cout << "Fun = " << Fun << "\n";
		}
//		std::cout << "Xup - Xdown = " << Xup - Xdown << "\n";		
//		std::cout << "DeltaX = " << DeltaX << "\n";
//		std::cout << "Integral = " << Integral << "\n";

		return Integral;
	}

	int ID_reserch(	dealii::Vector<dbl> &U_z, dealii::Vector<dbl> &x_coordinate_of_vertex,
					dealii::Vector<dbl> &y_coordinate_of_vertex, double Xcentr,
					double Ycentr, double rx, double ry, dealii::Point<2, double> point)						//Поиск вершину, близкой к центру (Xcentr,Ycentr)
	{
		double ID = 0;
		for(int i=0; i<U_z.size(); ++i)						//Пробежим по всем точкам, измеряя расстояние между данной точкой и центром
		{
			point(0) = x_coordinate_of_vertex(i);
			point(1) = y_coordinate_of_vertex(i);
//								if ( ( abs(point(0)-Xcentr) <= rx )&&( abs(point(1)-Ycentr) <= ry ) ) {ID = i;}
			if (   ( point(0)-Xcentr )*( point(0)-Xcentr ) + ( point(1)-Ycentr )*( point(1)-Ycentr ) <= rx*ry   ) {ID = i;}
		}
		return ID;
	}

	void Errors(	const int &i, dealii::Vector<dbl> &Vec1, dealii::Vector<dbl> &Vec2,
					double &deltABS, double &deltOTN,
					double &deltMaxABS, double &deltMaxOTN,
					const double &Xdown, const double &Xup,
					const double &Yleft, const double &Yright,
					const dealii::Point<2, double> &point, const double &h	)
	{
		if( (point(0)>Xdown)&&(point(0)<Xup) && (point(1)>Yleft)&&(point(1)<Yright)
				 && (point(0)!=h) && (point(0)!=2*h) )
		{
			deltABS = abs(Vec1(i)-Vec2(i));
			if(Vec2(i) != 0.0) {deltOTN = abs(deltABS / Vec2(i) * 100.0);}	else deltOTN = 0.0;
			if (deltMaxABS < deltABS) {deltMaxABS = deltABS;}
			if (deltMaxOTN < deltOTN) {deltMaxOTN = deltOTN;}
		}
/*
		deltABSU_x = abs(AU_x(i)-U_x(i));
		if(U_x(i) != 0.0) {deltOTNU_x = abs(deltABSU_x / U_x(i) * 100.0);}	else deltOTNU_x = 0.0;
		if (deltMaxABSU_x < deltABSU_x) {deltMaxABSU_x = deltABSU_x;}
		if (deltMaxOTNU_x < deltOTNU_x) {deltMaxOTNU_x = deltOTNU_x;}
*/
	}




	int main(		const double &Yung1,
					const double &Yung2,
					const double &Yung3,
					const double &Puasson1,
					const double &Puasson2,
					const double &Puasson3,
					const double &c0,
					const double &C0,
					const double &bb,
					const double &mu1,
					const double &mu2,
					const double &mu3,
					const double &h,
					const string FileOfGrid,
					const double &Yleft,
					const double &Yright,
					const double &Xup,
					const double &Xdown,
					const int &RefineNum,
					const double &IntegralCoefficient,
					const double &eps,
					FILE *FILES[]
			)
	{
		std::cout << "\n\nВызов функции main() из файла new_struct_elastic.cpp" << "\n\n";
		enum {x, y, z};
		FILE *FU_x			=	FILES[0];
		FILE *FU_y			=	FILES[1];
		FILE *FU_z			=	FILES[2];
		FILE *FU_x_grad		=	FILES[3];
		FILE *FU_y_grad		=	FILES[4];
		FILE *FU_z_grad		=	FILES[5];
		FILE *Ftau_xx		=	FILES[6];
		FILE *Ftau_yy		=	FILES[7];
		FILE *Ftau_xy		=	FILES[8];
		FILE *Ftau_zz		=	FILES[9];
		FILE *Ftau_zx		=	FILES[10];
		FILE *Ftau_zy		=	FILES[11];
		FILE *FAU_x			=	FILES[12];
		FILE *FAU_y			=	FILES[13];
		FILE *FAU_z			=	FILES[14];
		FILE *FAU_x_grad	=	FILES[15];
		FILE *FAU_y_grad	=	FILES[16];
		FILE *FAU_z_grad	=	FILES[17];
		FILE *FAtau_xx		=	FILES[18];
		FILE *FAtau_yy		=	FILES[19];
		FILE *FAtau_xy		=	FILES[20];
		FILE *FAtau_zz		=	FILES[21];
		FILE *FAtau_zx		=	FILES[22];
		FILE *FAtau_zy		=	FILES[23];


		//NIKOLA_ELASSTIC_PROBLEM
		if (1)
		{
			Domain<2> domain;

		    if (1)
			{
				cout << "\tYung1 	 = " << Yung1 << "\n";						//output
				cout << "\tYung2 	 = " << Yung2 << "\n";						//output
				cout << "\tYung3 	 = " << Yung3 << "\n";						//output
				cout << "\tPuasson1  = " << Puasson1 << "\n";					//output
				cout << "\tPuasson2  = " << Puasson2 << "\n";					//output
				cout << "\tPuasson3  = " << Puasson3 << "\n";					//output
				cout << "\tc0 		 = " << c0 << "\n";							//output
				cout << "\tC0 		 = " << C0 << "\n";							//output
				cout << "\tbb 		 = " << bb << "\n";							//output
				cout << "\tmu1 		 = " << mu1 << "\n";						//output
				cout << "\tmu2		 = " << mu2 << "\n";						//output
				cout << "\tmu3 		 = " << mu3 << "\n";						//output
				cout << "\tYleft 	 = " << Yleft << "\n";						//output
				cout << "\tYright 	 = " << Yright << "\n";						//output
				cout << "\tXup 		 = " << Xup << "\n";						//output
				cout << "\tXdown 	 = " << Xdown << "\n";						//output
				cout << "\tRefineNum = " << RefineNum << "\n";					//output
				cout << "\teps		 = " << eps << "\n";						//output
				dealii::GridIn<2> gridin;
				gridin.attach_triangulation(domain.grid);							//INPUT
				std::ifstream f(FileOfGrid);
				std::cout << "\tFileOfGrid = " << FileOfGrid << "\n";			//output
				Assert (dim==2, ExcInternalError());
				gridin.read_msh(f);
				domain.grid.refine_global(RefineNum);										//INPUT
					std::cout << "\tNumber of active cells:       "
				        	<< domain.grid.n_active_cells()
				        	<< std::endl;
		    }

		    dealii::FESystem<2,2> fe 
		        (dealii::FE_Q<2,2>(1), 2);
		    domain.dof_init (fe);

		    SystemsLinearAlgebraicEquations slae;
		    ATools ::trivial_prepare_system_equations (slae, domain);

		    LaplacianVector<2> element_matrix (domain.dof_handler.get_fe());
		    element_matrix.C .resize (3);											//один материал
		    EPTools ::set_isotropic_elascity{yung : Yung1, puasson : Puasson1}(element_matrix.C[0]);
		    EPTools ::set_isotropic_elascity{yung : Yung2, puasson : Puasson2}(element_matrix.C[1]);
		    EPTools ::set_isotropic_elascity{yung : Yung3, puasson : Puasson3}(element_matrix.C[2]);

			dealii::Vector<dbl> U_z(slae.solution.size() / 2);						//Вводим новый вектор
			dealii::Vector<dbl> tau_zx(slae.solution.size() / 2);						//Вводим новый вектор
			dealii::Vector<dbl> tau_zy(slae.solution.size() / 2);						//Вводим новый вектор
//			READ_FILE_TO_ARRAY::read_vector("out/FU_z_4225.gpl", U_z, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести
//			READ_FILE_TO_ARRAY::read_vector("out/Ftau_zx_4225.gpl", tau_zx, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести
//			READ_FILE_TO_ARRAY::read_vector("out/Ftau_zy_4225.gpl", tau_zy, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести

			READ_FILE_TO_ARRAY::read_vector("out/FU_z.gpl", U_z, 3, 4);				//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести
			READ_FILE_TO_ARRAY::read_vector("out/Ftau_zx.gpl", tau_zx, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести
			READ_FILE_TO_ARRAY::read_vector("out/Ftau_zy.gpl", tau_zy, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести

		    vec<arr<typename Nikola::SourceVector<2>::Func, 2>> U(3);
/*		    U[0][x] = [&U_z, &domain, &element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[0][x][x][z][z] * HCPTools::print_function_elastic<2>(U_z, domain, p);}; //Uz					//INPUT
		    U[0][y] = [&U_z, &domain, &element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[0][x][x][z][z] * HCPTools::print_function_elastic<2>(U_z, domain, p);}; //Uz					//INPUT
		    U[1][x] = [&U_z, &domain, &element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[1][x][x][z][z] * HCPTools::print_function_elastic<2>(U_z, domain, p);}; //Uz					//INPUT
		    U[1][y] = [&U_z, &domain, &element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[1][x][x][z][z] * HCPTools::print_function_elastic<2>(U_z, domain, p);}; //Uz					//INPUT
		    U[2][x] = [&U_z, &domain, &element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[2][x][x][z][z] * HCPTools::print_function_elastic<2>(U_z, domain, p);}; //Uz					//INPUT
		    U[2][y] = [&U_z, &domain, &element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[2][x][x][z][z] * HCPTools::print_function_elastic<2>(U_z, domain, p);}; //Uz					//INPUT
*/
		    U[0][x] = [&element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[0][x][x][z][z] * (-1.0)*(p(0)-c0);}; //Uz					//INPUT
		    U[0][y] = [&element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[0][x][x][z][z] * (-1.0)*(p(0)-c0);}; //Uz					//INPUT
		    U[1][x] = [&element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[1][x][x][z][z] * (-1.0)*(p(0)-c0);}; //Uz					//INPUT
		    U[1][y] = [&element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[1][x][x][z][z] * (-1.0)*(p(0)-c0);}; //Uz					//INPUT
		    U[2][x] = [&element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[2][x][x][z][z] * (-1.0)*(p(0)-c0);}; //Uz					//INPUT
		    U[2][y] = [&element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[2][x][x][z][z] * (-1.0)*(p(0)-c0);}; //Uz					//INPUT


		    vec<arr<typename Nikola::SourceVector<2>::Func, 2>> tau(3);
/*		    tau[0][x] = [&tau_zx, &domain] (const dealii::Point<2> &p) {return HCPTools::print_function_elastic<2>(tau_zx, domain, p);};				//tau_zx
		    tau[0][y] = [&tau_zy, &domain] (const dealii::Point<2> &p) {return HCPTools::print_function_elastic<2>(tau_zy, domain, p);};				//tau_zy
		    tau[1][x] = [&tau_zx, &domain] (const dealii::Point<2> &p) {return HCPTools::print_function_elastic<2>(tau_zx, domain, p);};				//tau_zx
		    tau[1][y] = [&tau_zy, &domain] (const dealii::Point<2> &p) {return HCPTools::print_function_elastic<2>(tau_zy, domain, p);};				//tau_zy
		    tau[2][x] = [&tau_zx, &domain] (const dealii::Point<2> &p) {return HCPTools::print_function_elastic<2>(tau_zx, domain, p);};				//tau_zx
		    tau[2][y] = [&tau_zy, &domain] (const dealii::Point<2> &p) {return HCPTools::print_function_elastic<2>(tau_zy, domain, p);};				//tau_zy
*/
		    tau[0][x] = [] (const dealii::Point<2> &p) {return 0.0;};	//tau_zx													//INPUT
		    tau[0][y] = [] (const dealii::Point<2> &p) {return 0.0;};	//tau_zy													//INPUT
		    tau[1][x] = [] (const dealii::Point<2> &p) {return 0.0;};	//tau_zx													//INPUT
		    tau[1][y] = [] (const dealii::Point<2> &p) {return 0.0;};	//tau_zy													//INPUT
		    tau[2][x] = [] (const dealii::Point<2> &p) {return 0.0;};	//tau_zx													//INPUT
		    tau[2][y] = [] (const dealii::Point<2> &p) {return 0.0;};	//tau_zy													//INPUT

		    Nikola::SourceVector<2> element_rhsv (U, tau, domain.dof_handler.get_fe());

		    Assembler ::assemble_matrix<2> (slae.matrix, element_matrix, domain.dof_handler);
		    Assembler ::assemble_rhsv<2> (slae.rhsv, element_rhsv, domain.dof_handler);


		    dealii::SolverControl solver_control (10000, eps);
		    dealii::SolverCG<> solver (solver_control);
		    solver.solve (
		            slae.matrix,
		            slae.solution,
		            slae.rhsv,
		            dealii::PreconditionIdentity()
		            );
			cout << "\n\n";





	//==============================================================================



			int quadrature_points = 0;												//Общее количество квадратурных точек
			int ID = -100;															//id вершины

		    dealii::Vector<dbl> x_coordinate_of_vertex(slae.solution.size() / 2);		//массив, содержит координаты x для каждой вершины сетки * 2 (т.к. 2 формулы U_x и U_y)
		    dealii::Vector<dbl> y_coordinate_of_vertex(slae.solution.size() / 2);		//массив, содержит координаты y для каждой вершины сетки * 2 (т.к. 2 формулы U_x и U_y)

		    dealii::Vector<dbl> U_x(slae.solution.size() / 2);						//хранит значения U_x
		    dealii::Vector<dbl> U_y(slae.solution.size() / 2);						//хранит значения U_y
		    dealii::Vector<dbl> U_x_gradX(slae.solution.size() / 2);				//хранит значения производных по x
		    dealii::Vector<dbl> U_y_gradX(slae.solution.size() / 2);				//хранит значения производных по x
		    dealii::Vector<dbl> U_x_gradY(slae.solution.size() / 2);				//хранит значения производных по x
		    dealii::Vector<dbl> U_y_gradY(slae.solution.size() / 2);				//хранит значения производных по x
		    dealii::Vector<dbl> tau_xx(slae.solution.size() / 2);						//массив, содержит значения tau_xx
		    dealii::Vector<dbl> tau_yy(slae.solution.size() / 2);						//массив, содержит значения tau_yy
		    dealii::Vector<dbl> tau_xy(slae.solution.size() / 2);						//массив, содержит значения tau_xy
		    dealii::Vector<dbl> tau_zz(slae.solution.size() / 2);						//массив, содержит значения tau_xy

							dealii::Vector<dbl> AU_x(slae.solution.size() / 2);						//хранит значения AU_x
							dealii::Vector<dbl> AU_y(slae.solution.size() / 2);						//хранит значения AU_y
							dealii::Vector<dbl> AU_x_gradX(slae.solution.size() / 2);				//хранит значения производных по x
							dealii::Vector<dbl> AU_y_gradX(slae.solution.size() / 2);				//хранит значения производных по x
							dealii::Vector<dbl> AU_x_gradY(slae.solution.size() / 2);				//хранит значения производных по x
							dealii::Vector<dbl> AU_y_gradY(slae.solution.size() / 2);				//хранит значения производных по x
							dealii::Vector<dbl> Atau_xx(slae.solution.size() / 2);						//массив, содержит значения Atau_xx
							dealii::Vector<dbl> Atau_yy(slae.solution.size() / 2);						//массив, содержит значения Atau_yy
							dealii::Vector<dbl> Atau_xy(slae.solution.size() / 2);						//массив, содержит значения Atau_xy
							dealii::Vector<dbl> Atau_zz(slae.solution.size() / 2);						//массив, содержит значения Atau_xy



		    {
		        dealii::QGauss<2>  quadrature_formula(2);

		        dealii::FEValues<2> fe_values (domain.dof_handler.get_fe(), quadrature_formula,
		                dealii::update_quadrature_points | dealii::update_gradients | dealii::update_JxW_values |
		                dealii::update_values);

		        cst dofs_per_cell = fe.dofs_per_cell;								// 8
		        cst n_q_points = quadrature_formula.size();							// 4, возвращает кол-во квадратурных точек


		        vec<dealii::types::global_dof_index> local_dof_indices (dofs_per_cell);

				int number_of_dots = 0;
		        auto cell = domain.dof_handler.begin_active();
		        auto endc = domain.dof_handler.end();
				//Пробежать по всем ячейкам, узнать номер DoF, 
				cell = domain.dof_handler.begin_active();
			    endc = domain.dof_handler.end();



				

	//------------------------------------------------------------------------------
	//------------------------------------------------------------------------------
	//Проверка работы функции градиента
	/*
			    for (; cell != endc; ++cell)
	   	        {
	//cell->vertex_dof_index (vertex, i)		-		возвращает глобальный DoF-индекс для вершины vertex,
	//													степени i. i - для выбора функции (Например функции, 0 - U_x, 1 - U_y)
	//ID = cell->vertex_dof_index (0, 0);		-		Присваивание буферной переменной глобального номера DoF для первой функции
	//test_length(ID) = ID;		-		Присваивание ID-ому элементу массива значения своего собственного номера ID
	//ID = cell->vertex_dof_index (0, 1);		-		Присваивание буферной переменной глобального номера DoF для второй функции
	//test_length(ID) = ID;		-		Присваивание ID-ому элементу массива значения своего собственного номера ID
	//x_coordinate_of_vertex(ID) = cell->vertex(0)[0];		-		Присваивание ID-ому элементу массива значения x-координат вершины
	//y_coordinate_of_vertex(ID) = cell->vertex(0)[1];		-		Присваивание ID-ому элементу массива значения y-координат вершины
					ID = cell->vertex_dof_index (0, 0);
					test_length(ID) = ID;
					x_coordinate_of_vertex(ID) = cell->vertex(0)[0];
					y_coordinate_of_vertex(ID) = cell->vertex(0)[1];
					ID = cell->vertex_dof_index (0, 1);
					test_length(ID) = ID;
					x_coordinate_of_vertex(ID) = cell->vertex(0)[0];
					y_coordinate_of_vertex(ID) = cell->vertex(0)[1];
						ID = cell->vertex_dof_index (1, 0);
						test_length(ID) = ID;
						x_coordinate_of_vertex(ID) = cell->vertex(1)[0];
						y_coordinate_of_vertex(ID) = cell->vertex(1)[1];
						ID = cell->vertex_dof_index (1, 1);
						test_length(ID) = ID;
						x_coordinate_of_vertex(ID) = cell->vertex(1)[0];
						y_coordinate_of_vertex(ID) = cell->vertex(1)[1];
							ID = cell->vertex_dof_index (2, 0);
							test_length(ID) = ID;
							x_coordinate_of_vertex(ID) = cell->vertex(2)[0];
							y_coordinate_of_vertex(ID) = cell->vertex(2)[1];
							ID = cell->vertex_dof_index (2, 1);
							test_length(ID) = ID;
							x_coordinate_of_vertex(ID) = cell->vertex(2)[0];
							y_coordinate_of_vertex(ID) = cell->vertex(2)[1];
								ID = cell->vertex_dof_index (3, 0);
								test_length(ID) = ID;
								x_coordinate_of_vertex(ID) = cell->vertex(3)[0];
								y_coordinate_of_vertex(ID) = cell->vertex(3)[1];
								ID = cell->vertex_dof_index (3, 1);
								test_length(ID) = ID;
								x_coordinate_of_vertex(ID) = cell->vertex(3)[0];
								y_coordinate_of_vertex(ID) = cell->vertex(3)[1];
				}
				printf("\t\tnumber_of_dots = %8d\n", number_of_dots);


	//создание параболы для проверки работы функции градиента
				std::cout << "domain.dof_handler.n_dofs() = " << domain.dof_handler.n_dofs() << "\n";
				for(int i=0; i<domain.dof_handler.n_dofs(); ++i)
				{
	//				U_x(i) = x_coordinate_of_vertex(i) * x_coordinate_of_vertex(i) + y_coordinate_of_vertex(i) * y_coordinate_of_vertex(i);
					slae.solution(i) = sin(x_coordinate_of_vertex(i)) + cos(y_coordinate_of_vertex(i));
	//				slae.solution(i) = x_coordinate_of_vertex(i) * x_coordinate_of_vertex(i) + y_coordinate_of_vertex(i) * y_coordinate_of_vertex(i);
	//				slae.solution(i) = x_coordinate_of_vertex(i)  + y_coordinate_of_vertex(i);
					std::cout << "   x = " << x_coordinate_of_vertex(i) << "\t\t";
					std::cout << "   y = " << y_coordinate_of_vertex(i) << "\t\t";
					std::cout << "slae.solution(" << i << ") = " << slae.solution(i) << "\n";
				}

	//Вывод в файл
				FILE *TestFile1;
				TestFile1 = fopen ("out/TestFile1.gpl", "w");
				for(int i=0; i<domain.dof_handler.n_dofs(); ++i)
				{
	//																					№			"x"						"y"				  			"F"						"dX"					"dY"
					fprintf(  TestFile1, "%8d\t%5.35f\t%5.35f\t%5.35f\t%5.35f\t%5.35f\r\n",	i, x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), slae.solution(i), cos(x_coordinate_of_vertex(i)), -sin(y_coordinate_of_vertex(i)) );
	//				fprintf(  TestFile1, "%8d\t%5.35f\t%5.35f\t%5.35f\t%5.35f\t%5.35f\r\n",	i, x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), slae.solution(i), 2*x_coordinate_of_vertex(i), 2*y_coordinate_of_vertex(i) );
	//				fprintf(  TestFile1, "%8d\t%5.35f\t%5.35f\t%5.35f\t%5.35f\t%5.35f\r\n",	i, x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), slae.solution(i), 				1.0, 					1.0 );
				}
				fclose(TestFile1);
	*/
	//------------------------------------------------------------------------------
	//------------------------------------------------------------------------------

	//------------------------------------------------------------------------------
	//инедексация векторов 
				int D = 0.0;
			    for (; cell != endc; ++cell)
	   	        {
	//cell->vertex_dof_index (vertex, i)		-		возвращает глобальный DoF-индекс для вершины vertex,
	//													степени i. i - для выбора функции (Например функции, 0 - U_x, 1 - U_y)
	//ID = cell->vertex_dof_index (0, 0);		-		Присваивание буферной переменной глобального номера DoF для первой функции
	//test_length(ID) = ID;		-		Присваивание ID-ому элементу массива значения своего собственного номера ID
	//ID = cell->vertex_dof_index (0, 1);		-		Присваивание буферной переменной глобального номера DoF для второй функции
	//test_length(ID) = ID;		-		Присваивание ID-ому элементу массива значения своего собственного номера ID
	//x_coordinate_of_vertex(ID) = cell->vertex(0)[0];		-		Присваивание ID-ому элементу массива значения x-координат вершины
	//y_coordinate_of_vertex(ID) = cell->vertex(0)[1];		-		Присваивание ID-ому элементу массива значения y-координат вершины
					ID = cell->vertex_dof_index (0, 0);
					D = ID / 2;
					x_coordinate_of_vertex(D) = cell->vertex(0)[0];
					y_coordinate_of_vertex(D) = cell->vertex(0)[1];
					ID = cell->vertex_dof_index (0, 1);
					D = ID / 2;
					x_coordinate_of_vertex(D) = cell->vertex(0)[0];
					y_coordinate_of_vertex(D) = cell->vertex(0)[1];
						ID = cell->vertex_dof_index (1, 0);
						D = ID / 2;
						x_coordinate_of_vertex(D) = cell->vertex(1)[0];
						y_coordinate_of_vertex(D) = cell->vertex(1)[1];
						ID = cell->vertex_dof_index (1, 1);
						D = ID / 2;
						x_coordinate_of_vertex(D) = cell->vertex(1)[0];
						y_coordinate_of_vertex(D) = cell->vertex(1)[1];
							ID = cell->vertex_dof_index (2, 0);
							D = ID / 2;
							x_coordinate_of_vertex(D) = cell->vertex(2)[0];
							y_coordinate_of_vertex(D) = cell->vertex(2)[1];
							ID = cell->vertex_dof_index (2, 1);
							D = ID / 2;
							x_coordinate_of_vertex(D) = cell->vertex(2)[0];
							y_coordinate_of_vertex(D) = cell->vertex(2)[1];
								ID = cell->vertex_dof_index (3, 0);
								D = ID / 2;
								x_coordinate_of_vertex(D) = cell->vertex(3)[0];
								y_coordinate_of_vertex(D) = cell->vertex(3)[1];
								ID = cell->vertex_dof_index (3, 1);
								D = ID / 2;
								x_coordinate_of_vertex(D) = cell->vertex(3)[0];
								y_coordinate_of_vertex(D) = cell->vertex(3)[1];
				}
//				printf("\t\tnumber_of_dots = %8d\n", number_of_dots);

	//------------------------------------------------------------------------------


/*
//------------------------------------------------------------------------------
				for(int i=0; i<slae.solution.size(); ++i)
				{
					dealii::Point<2, double> point(2.0,444.4);
					point(0) = x_coordinate_of_vertex(i / 2);
					point(1) = y_coordinate_of_vertex(i / 2);
					slae.solution(i) = -point(0);
				}
//------------------------------------------------------------------------------
*/





	//------------------------------------------------------------------------------
	//Создание файла GRAD.gpl и 
	//создать файл GRAD.gpl - содержит производные
	//			HCPTools::print_heat_gradient<2>(slae.solution, domain, "GRAD");
	//			std::cout << "HCPTools::print_heat_gradient<2> \n";
				HCPTools::print_elasticity_gradient<2>(slae.solution, domain, "GRAD");			//Вызов градиент-функции для U_x
				std::cout << "HCPTools::print_elasticity_gradient<2> \n";

				int w = 0;																	//w - для результата работы функции read_files()
				int M = slae.solution.size(), N = 6;										//M - количество строк в файле, N - количество столбцов в файле. M = 2178, N = 5;
				w = READ_FILE_TO_ARRAY::read_files_for_elelastic(	x_coordinate_of_vertex, y_coordinate_of_vertex,
													U_x, U_y,
													U_x_gradX, U_x_gradY,
													U_y_gradX, U_y_gradY, N);

				std::cout << "slae.solution.size() = " << slae.solution.size() << "\n";
				std::cout << "U_x.size() = " << U_x.size() << "\n";
				std::cout << "U_y.size() = " << U_y.size() << "\n";
				std::cout << "U_x_gradX.size() = " << U_x_gradX.size() << "\n";
				std::cout << "U_x_gradY.size() = " << U_x_gradY.size() << "\n";
				std::cout << "U_y_gradX.size() = " << U_y_gradX.size() << "\n";
				std::cout << "U_y_gradY.size() = " << U_y_gradY.size() << "\n";
				std::cout << "\n";
	//------------------------------------------------------------------------------
			    
















	//==============================================================================
	//==============================================================================

				for(int i=0; i<slae.solution.size() / 2; ++i)
				{
	//Для того чтобы подставлять значения кординат точек
					dealii::Point<2, double> point(2.0,444.4);
					point(0) = x_coordinate_of_vertex(i);
					point(1) = y_coordinate_of_vertex(i);

	//-------------------------------------			//INPUT

					if(point(0) < h)
					{
						AU_x(i) = Puasson1 / 2 * ( -point(1)*point(1) + (point(0) - c0)*(point(0) - c0) + 0*2*C0 + (bb*bb-1)/12 );
					}
					else if(point(0) < 2*h)
					{
						AU_x(i) = Puasson2 / 2 * ( -point(1)*point(1) + (point(0) - c0)*(point(0) - c0) + 0*2*C0 + (bb*bb-1)/12 );
					}
					else
					{
						AU_x(i) = Puasson3 / 2 * ( -point(1)*point(1) + (point(0) - c0)*(point(0) - c0) + 0*2*C0 + (bb*bb-1)/12 );
					}
					fprintf( FAU_x, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), AU_x(i) );

				}
	//==============================================================================
	//==============================================================================


	//------------------------------------------------------------------------------
	//Вычисление интеграла
			{
				double AddConstant = 0.0;
				double integralConst = 0.0;
				cout << " => (in beginning) Integral(AU_x) = " << Integral(AU_x, domain, fe_values, n_q_points, dofs_per_cell) << "\n";

				for(int j = 0; j < 3; j++)
				{
					integralConst = 0.95 * Integral(AU_x, domain, fe_values, n_q_points, dofs_per_cell) / bb;
					Add_to_vector(AU_x, -integralConst);
					AddConstant += -integralConst;
				}
				integralConst = Integral(AU_x, domain, fe_values, n_q_points, dofs_per_cell);
				cout << " => Integral(AU_x) = " << integralConst << "\n";
				cout << "AddConstant = " << AddConstant << "\n";
				std::cout << "\n";
				std::cout << "\n";
			}
	//------------------------------------------------------------------------------











	//==============================================================================
	//==============================================================================
	//цикл вывода аналитических функций по всем вершинам
				for(int i=0; i<slae.solution.size() / 2; ++i)
				{
	//Для того чтобы подставлять значения кординат точек
					dealii::Point<2, double> point(2.0,444.4);
					point(0) = x_coordinate_of_vertex(i);
					point(1) = y_coordinate_of_vertex(i);

	//-------------------------------------			//INPUT
/*
					if(point(0) < h)
					{
						AU_x(i) = Puasson1 / 2 * ( -point(1)*point(1) + (point(0) - c0)*(point(0) - c0) + 0*2*C0 + (bb*bb*bb-bb)/12 );
					}
					else if(point(0) < 2*h)
					{
						AU_x(i) = Puasson2 / 2 * ( -point(1)*point(1) + (point(0) - c0)*(point(0) - c0) + 0*2*C0 + (bb*bb*bb-bb)/12 );
					}
					else
					{
						AU_x(i) = Puasson3 / 2 * ( -point(1)*point(1) + (point(0) - c0)*(point(0) - c0) + 0*2*C0 + (bb*bb*bb-bb)/12 );
					}
					fprintf( FAU_x, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), AU_x(i) );
*/

					if(point(0) < h)
					{
						AU_y(i) = Puasson1 * (point(0) - c0) * point(1);
					}
					else if(point(0) < 2*h)
					{
						AU_y(i) = Puasson2 * (point(0) - c0) * point(1);
					}
					else
					{
						AU_y(i) = Puasson3 * (point(0) - c0) * point(1);
					}
					fprintf( FAU_y, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), AU_y(i) );

	//-------------------------------------			//INPUT

					if(point(0) < h)
					{
						AU_x_gradX(i) = Puasson1 / 2 * 2 * (point(0) - c0);
						AU_x_gradY(i) = Puasson1 / 2 * (-2) * point(1);
					}
					else if(point(0) < 2*h)
					{
						AU_x_gradX(i) = Puasson2 / 2 * 2 * (point(0) - c0);
						AU_x_gradY(i) = Puasson2 / 2 * (-2) * point(1);
					}
					else
					{
						AU_x_gradX(i) = Puasson3 / 2 * 2 * (point(0) - c0);
						AU_x_gradY(i) = Puasson3 / 2 * (-2) * point(1);
					}
					fprintf( FAU_x_grad, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1),
																AU_x(i), AU_x_gradX(i), AU_x_gradY(i) );

					if(point(0) < h)
					{
						AU_y_gradX(i) = Puasson1 * point(1);
						AU_y_gradY(i) = Puasson1 * (point(0) - c0);
					}
					else if(point(0) < 2*h)
					{
						AU_y_gradX(i) = Puasson2 * point(1);
						AU_y_gradY(i) = Puasson2 * (point(0) - c0);
					}
					else
					{
						AU_y_gradX(i) = Puasson3 * point(1);
						AU_y_gradY(i) = Puasson3 * (point(0) - c0);
					}
					fprintf( FAU_y_grad, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1),
																AU_y(i), AU_y_gradX(i), AU_y_gradY(i) );

	//-------------------------------------

					Atau_xx(i) = 0.0;
					fprintf( FAtau_xx, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), Atau_xx(i) );

					Atau_yy(i) = 0.0;
					fprintf( FAtau_yy, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), Atau_yy(i) );

					Atau_xy(i) = 0.0;
					fprintf( FAtau_xy, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), Atau_xy(i) );

	//-------------------------------------			//INPUT

					if(point(0) < h)
					{
						Atau_zz(i) = -Yung1*(point(0) - c0);
					}
					else if(point(0) < 2*h)
					{
						Atau_zz(i) = -Yung2*(point(0) - c0);
					}
					else
					{
						Atau_zz(i) = -Yung3*(point(0) - c0);
					}
					fprintf( FAtau_zz, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), Atau_zz(i) );

				}

	//==============================================================================
	//==============================================================================















/*
	//------------------------------------------------------------------------------
	//Вычисление интеграла
			{
				//находим примерный радиус пятна для поиска центральной точки
				double AddConstant = 0.0;
				double integralConst = 0.0;
				double ID = 0;
				double AveArea = (Xup - Xdown)*(Yright - Yleft) / U_x.size(); 				//Найдём среднюю площадь на одну вершину. Для этого общуую площадь разделим на кол-во вершин.
				double rx = sqrt(AveArea) / 50;
				double ry = sqrt(AveArea) / 50;
				double Xcentr = (Xup + Xdown) / 2;											//Найдём центр 
				double Ycentr = (Yright + Yleft) / 2;										//Найдём центр 
				dealii::Point<2, double> point(2.0,444.4);
				//ищем № центральной точки пробегая по всем вершинам, если не находим - увичеваем радиус пятна в 1.1 раз
						cout << "AveArea = " << AveArea << "\n";																		//output
						cout << "rx = " << rx << "\n";																					//output
						cout << "ry = " << ry << "\n";																					//output
						ID = ID_reserch(U_x, x_coordinate_of_vertex, y_coordinate_of_vertex, Xcentr, Ycentr, rx, ry, point);
						while(ID == 0)
						{
							iter++;
							rx *= 1.1; ry *= 1.1;
							ID = ID_reserch(U_x, x_coordinate_of_vertex, y_coordinate_of_vertex, Xcentr, Ycentr, rx, ry, point);
						}
						cout << " x_coordinate_of_vertex(ID) = " <<  x_coordinate_of_vertex(ID) << "\n";								//output
						cout << " y_coordinate_of_vertex(ID) = " <<  y_coordinate_of_vertex(ID) << "\n";								//output
						cout << " Xcentr = " <<  Xcentr << "\n";																		//output
						cout << " Ycentr = " <<  Ycentr << "\n";																		//output
						cout << "iter = " << iter << "\n";																				//output
						cout << "ID = " << ID << "\n";																					//output
						cout << "U_x(ID) = " << U_x(ID) << "\n";
						cout << "AU_x(ID) = " << AU_x(ID) << "\n";
						cout << "U_x(ID) - AU_x(ID) = " << U_x(ID) - AU_x(ID) << "\n";
						std::cout << "\n";
				//Высчитываем разницу значений функций в точке №(ID). Приближаем функцию AU_x к функции U_x с каждым шагом.
								for(int j = 0; j < 30; j++)
								{
									integralConst =  0.5*( U_x(ID) - AU_x(ID) );
									Add_to_vector(U_x, -integralConst);
									AddConstant += -integralConst;
									cout << "U_x(ID) - AU_x(ID) = " << U_x(ID) - AU_x(ID) << "\n";
								}
								integralConst = Integral(U_x, domain, fe_values, n_q_points, dofs_per_cell);	
								cout << " => Integral(U_x) = " << integralConst << "\n";
								cout << "AddConstant = " << AddConstant << "\n";

								cout << " => Integral(AU_x) = " << Integral(AU_x, domain, fe_values, n_q_points, dofs_per_cell) << "\n";
								std::cout << "\n";

			}
	//------------------------------------------------------------------------------
*/




	//------------------------------------------------------------------------------
	//Вычисление интеграла
			{
				double AddConstant = 0.0;
				double integralConst = 0.0;
				cout << " => (in beginning) Integral(U_x) = " << Integral(U_x, domain, fe_values, n_q_points, dofs_per_cell) << "\n";

				for(int j = 0; j < 3; j++)
				{
					integralConst = 0.95 * Integral(U_x, domain, fe_values, n_q_points, dofs_per_cell) / bb;
					Add_to_vector(U_x, -integralConst);
					AddConstant += -integralConst;
				}
				integralConst = Integral(U_x, domain, fe_values, n_q_points, dofs_per_cell);	
				cout << " => Integral(U_x) = " << integralConst << "\n";
				cout << "AddConstant = " << AddConstant << "\n";
				std::cout << "\n";
				std::cout << "\n";
			}
	//------------------------------------------------------------------------------













	//==============================================================================
	//==============================================================================

	double		deltABSU_x = -100;
	double		deltOTNU_x = -100;
	double		deltMaxABSU_x = -100;
	double		deltMaxOTNU_x = -100;

	double		deltABSU_x_gradX = -100;
	double		deltOTNU_x_gradX = -100;
	double		deltMaxABSU_x_gradX = -100;
	double		deltMaxOTNU_x_gradX = -100;

	double		deltABSU_x_gradY = -100;
	double		deltOTNU_x_gradY = -100;
	double		deltMaxABSU_x_gradY = -100;
	double		deltMaxOTNU_x_gradY = -100;

	double		deltABSU_y = -100;
	double		deltOTNU_y = -100;
	double		deltMaxABSU_y = -100;
	double		deltMaxOTNU_y = -100;

	double		deltABSU_y_gradX = -100;
	double		deltOTNU_y_gradX = -100;
	double		deltMaxABSU_y_gradX = -100;
	double		deltMaxOTNU_y_gradX = -100;

	double		deltABSU_y_gradY = -100;
	double		deltOTNU_y_gradY = -100;
	double		deltMaxABSU_y_gradY = -100;
	double		deltMaxOTNU_y_gradY = -100;

	double		deltABStau_xx = -100;
	double		deltAVEtau_xx = 0.0;
	double		deltOTNtau_xx = -100;
	double		deltMaxABStau_xx = -100;
	double		deltMaxOTNtau_xx = -100;

	double		deltABStau_yy = -100;
	double		deltAVEtau_yy = 0.0;
	double		deltOTNtau_yy = -100;
	double		deltMaxABStau_yy = -100;
	double		deltMaxOTNtau_yy = -100;

	double		deltABStau_xy = -100;
	double		deltAVEtau_xy = 0.0;
	double		deltOTNtau_xy = -100;
	double		deltMaxABStau_xy = -100;
	double		deltMaxOTNtau_xy = -100;

	double		deltABStau_zz = -100;
	double		deltAVEtau_zz = 0.0;
	double		deltOTNtau_zz = -100;
	double		deltMaxABStau_zz = -100;
	double		deltMaxOTNtau_zz = -100;

	//тот самых цикл вывода по всем вершинам
				for(int i=0; i<slae.solution.size() / 2; ++i)
				{
	//Для того чтобы подставлять значения кординат точек
					dealii::Point<2, double> point(2.0,444.4);
					point(0) = x_coordinate_of_vertex(i);
					point(1) = y_coordinate_of_vertex(i);

	//-------------------------------------

					fprintf( FU_x, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), U_x(i) );
					fprintf( FU_y, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), U_y(i) );

	//-------------------------------------

					fprintf( FU_x_grad, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), U_x(i), U_x_gradX(i), U_x_gradY(i) );
					Errors(	i, AU_x, U_x, deltABSU_x, deltOTNU_x, deltMaxABSU_x, deltMaxOTNU_x, Xdown, Xup, Yleft, Yright, point, h	);
					Errors(	i, AU_x_gradX, U_x_gradX, deltABSU_x_gradX, deltOTNU_x_gradX, deltMaxABSU_x_gradX, deltMaxOTNU_x_gradX, Xdown, Xup, Yleft, Yright, point, h	);
					Errors(	i, AU_x_gradY, U_x_gradY, deltABSU_x_gradY, deltOTNU_x_gradY, deltMaxABSU_x_gradY, deltMaxOTNU_x_gradY, Xdown, Xup, Yleft, Yright, point, h	);

					fprintf( FU_y_grad, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), U_y(i), U_y_gradX(i), U_y_gradY(i) );
					Errors(	i, AU_y, U_y, deltABSU_y, deltOTNU_y, deltMaxABSU_y, deltMaxOTNU_y, Xdown, Xup, Yleft, Yright, point, h	);
					Errors(	i, AU_y_gradX, U_y_gradX, deltABSU_y_gradX, deltOTNU_y_gradX, deltMaxABSU_y_gradX, deltMaxOTNU_y_gradX, Xdown, Xup, Yleft, Yright, point, h	);
					Errors(	i, AU_y_gradY, U_y_gradY, deltABSU_y_gradY, deltOTNU_y_gradY, deltMaxABSU_y_gradY, deltMaxOTNU_y_gradY, Xdown, Xup, Yleft, Yright, point, h	);

	//-------------------------------------			//INPUT

					if(point(0) < h)
					{
						tau_xx(i)	=	element_matrix.C[0][x][x][x][x] * U_x_gradX(i)
									+	element_matrix.C[0][x][x][y][y] * U_y_gradY(i)
									+	U[0][x](point);	//Разделим, т. к. при вызове функции возвращается значение,
	//								+	element_matrix.C[0][x][x][z][z] * U[0][x](point) / element_matrix.C[0][x][x][z][z];	//Разделим, т. к. при вызове функции возвращается значение,
																															//умноженное на element_matrix.C[0][x][x][z][z];

						tau_yy(i)	=	element_matrix.C[0][y][y][x][x] * U_x_gradX(i)
									+	element_matrix.C[0][y][y][y][y] * U_y_gradY(i)
									+	U[0][y](point);	//Разделим, т. к. при вызове функции возвращается значение,
	//								+	element_matrix.C[0][y][y][z][z] * U[0][y](point) / element_matrix.C[0][x][x][z][z];	//Разделим, т. к. при вызове функции возвращается значение,
																														//умноженное на element_matrix.C[0][x][x][z][z];
					}
					else if(point(0) < 2*h)
					{
						tau_xx(i)	=	element_matrix.C[1][x][x][x][x] * U_x_gradX(i)
									+	element_matrix.C[1][x][x][y][y] * U_y_gradY(i)
									+	U[1][x](point);	//Разделим, т. к. при вызове функции возвращается значение,
	//								+	element_matrix.C[1][x][x][z][z] * U[1][x](point) / element_matrix.C[1][x][x][z][z];	//Разделим, т. к. при вызове функции возвращается значение,
																															//умноженное на element_matrix.C[1][x][x][z][z];

						tau_yy(i)	=	element_matrix.C[1][y][y][x][x] * U_x_gradX(i)
									+	element_matrix.C[1][y][y][y][y] * U_y_gradY(i)
									+	U[1][y](point);	//Разделим, т. к. при вызове функции возвращается значение,
	//								+	element_matrix.C[1][y][y][z][z] * U[1][y](point) / element_matrix.C[1][x][x][z][z];	//Разделим, т. к. при вызове функции возвращается значение,
																														//умноженное на element_matrix.C[1][x][x][z][z];
					}
					else
					{
						tau_xx(i)	=	element_matrix.C[2][x][x][x][x] * U_x_gradX(i)
									+	element_matrix.C[2][x][x][y][y] * U_y_gradY(i)
									+	U[2][x](point);	//Разделим, т. к. при вызове функции возвращается значение,
	//								+	element_matrix.C[2][x][x][z][z] * U[2][x](point) / element_matrix.C[2][x][x][z][z];	//Разделим, т. к. при вызове функции возвращается значение,
																															//умноженное на element_matrix.C[2][x][x][z][z];

						tau_yy(i)	=	element_matrix.C[2][y][y][x][x] * U_x_gradX(i)
									+	element_matrix.C[2][y][y][y][y] * U_y_gradY(i)
									+	U[2][y](point);	//Разделим, т. к. при вызове функции возвращается значение,
	//								+	element_matrix.C[2][y][y][z][z] * U[1][y](point) / element_matrix.C[2][x][x][z][z];	//Разделим, т. к. при вызове функции возвращается значение,
																															//умноженное на element_matrix.C[2][x][x][z][z];
					}


					if(point(0) < h)
					{
						tau_xy(i)	=	element_matrix.C[0][x][y][x][y] * 
										( U_x_gradY(i) + U_y_gradX(i) );
					}
					else if(point(0) < 2*h)
					{
						tau_xy(i)	=	element_matrix.C[1][x][y][x][y] * 
										( U_x_gradY(i) + U_y_gradX(i) );
					}
					else
					{
						tau_xy(i)	=	element_matrix.C[2][x][y][x][y] * 
										( U_x_gradY(i) + U_y_gradX(i) );
					}
				

					fprintf( Ftau_xx, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), tau_xx(i) );
					fprintf( Ftau_yy, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), tau_yy(i) );
					fprintf( Ftau_xy, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), tau_xy(i) );
					Errors(	i, Atau_xx, tau_xx, deltABStau_xx, deltOTNtau_xx, deltMaxABStau_xx, deltMaxOTNtau_xx, Xdown, Xup, Yleft, Yright, point, h	);
					Errors(	i, Atau_yy, tau_yy, deltABStau_yy, deltOTNtau_yy, deltMaxABStau_yy, deltMaxOTNtau_yy, Xdown, Xup, Yleft, Yright, point, h	);
					Errors(	i, Atau_xy, tau_xy, deltABStau_xy, deltOTNtau_xy, deltMaxABStau_xy, deltMaxOTNtau_xy, Xdown, Xup, Yleft, Yright, point, h	);

	//-------------------------------------			//INPUT

					if(point(0) < h)
					{
						tau_zz(i)	=	element_matrix.C[0][x][z][x][z] * tau_xx(i)
									+	element_matrix.C[0][y][z][y][z] * tau_yy(i)
									+	Yung1 * U[0][x](point) / element_matrix.C[0][x][x][z][z];	//Разделим, т. к. при вызове функции возвращается значение,
																									//умноженное на element_matrix.C[0][x][x][z][z]
					}
					else if(point(0) < 2*h)
					{
						tau_zz(i)	=	element_matrix.C[1][x][z][x][z] * tau_xx(i)
									+	element_matrix.C[1][y][z][y][z] * tau_yy(i)
									+	Yung2 * U[1][x](point) / element_matrix.C[1][x][x][z][z];	//Разделим, т. к. при вызове функции возвращается значение,
																									//умноженное на element_matrix.C[1][x][x][z][z]
					}
					else
					{
						tau_zz(i)	=	element_matrix.C[2][x][z][x][z] * tau_xx(i)
									+	element_matrix.C[2][y][z][y][z] * tau_yy(i)
									+	Yung3 * U[2][x](point) / element_matrix.C[2][x][x][z][z];	//Разделим, т. к. при вызове функции возвращается значение,
																									//умноженное на element_matrix.C[2][x][x][z][z]
					}
					fprintf( Ftau_zz, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1),  tau_zz(i));
					Errors(	i, Atau_zz, tau_zz, deltABStau_zz, deltOTNtau_zz, deltMaxABStau_zz, deltMaxOTNtau_zz, Xdown, Xup, Yleft, Yright, point, h	);

				}

				std::cout << "deltMaxABSU_x			= " << deltMaxABSU_x << "\t\t";
				std::cout << "deltMaxOTNU_x			= " << deltMaxOTNU_x << "%\n";

				std::cout << "deltMaxABSU_x_gradX		= " << deltMaxABSU_x_gradX << "\t\t";
				std::cout << "deltMaxOTNU_x_gradX	= " << deltMaxOTNU_x_gradX << "%\n";

				std::cout << "deltMaxABSU_x_gradY		= " << deltMaxABSU_x_gradY << "\t\t";
				std::cout << "deltMaxOTNU_x_gradY	= " << deltMaxOTNU_x_gradY << "%\n";
				std::cout << "\n";


				std::cout << "deltMaxABSU_y			= " << deltMaxABSU_y << "\t\t";
				std::cout << "deltMaxOTNU_y			= " << deltMaxOTNU_y << "%\n";

				std::cout << "deltMaxABSU_y_gradX		= " << deltMaxABSU_y_gradX << "\t\t";
				std::cout << "deltMaxOTNU_y_gradX	= " << deltMaxOTNU_y_gradX << "%\n";

				std::cout << "deltMaxABSU_y_gradY		= " << deltMaxABSU_y_gradY << "\t\t";
				std::cout << "deltMaxOTNU_y_gradY	= " << deltMaxOTNU_y_gradY << "%\n";
				std::cout << "\n";


				std::cout << "deltMaxABStau_xx		= " << deltMaxABStau_xx << "\t\t\t";
				std::cout << "deltMaxOTNtau_xx		= " << deltMaxOTNtau_xx << "%\t\t\t";
				std::cout << "deltAVEtau_xx			= " << deltAVEtau_xx / (slae.solution.size() / 2) << "\n";

				std::cout << "deltMaxABStau_yy		= " << deltMaxABStau_yy << "\t\t\t";
				std::cout << "deltMaxOTNtau_yy		= " << deltMaxOTNtau_yy << "%\t\t\t";
				std::cout << "deltAVEtau_yy			= " << deltAVEtau_yy / (slae.solution.size() / 2) << "\n";

				std::cout << "deltMaxABStau_xy		= " << deltMaxABStau_xy << "\t\t\t";
				std::cout << "deltMaxOTNtau_xy		= " << deltMaxOTNtau_xy << "%\t\t\t";
				std::cout << "deltAVEtau_xy			= " << deltAVEtau_xy / (slae.solution.size() / 2) << "\n";

				std::cout << "deltMaxABStau_zz		= " << deltMaxABStau_zz << "\t\t\t";
				std::cout << "deltMaxOTNtau_zz		= " << deltMaxOTNtau_zz << "%\t\t\t";
				std::cout << "deltAVEtau_zz			= " << deltAVEtau_zz / (slae.solution.size() / 2) << "\n";

				std::cout << "\tFileOfGrid			= " << FileOfGrid << "_" << domain.grid.n_active_cells() << "_elastic" <<"\n";
	//==============================================================================
	//==============================================================================

				{

					int N = 500;
					double IntegralOfTau_xx = 0.0;
					IntegralOfTau_xx = IntegralOfCoordinateX(tau_xx, domain, fe_values, n_q_points, dofs_per_cell, 0.0, Yleft, Yright, N);
					cout << "IntegralOfTau_xx" << "(0.0) = " << IntegralOfTau_xx << "\n";

					IntegralOfTau_xx = IntegralOfCoordinateX(tau_xx, domain, fe_values, n_q_points, dofs_per_cell, 0.166, Yleft, Yright, N);
					cout << "IntegralOfTau_xx" << "(0.166) = " << IntegralOfTau_xx << "\n";

					IntegralOfTau_xx = IntegralOfCoordinateX(tau_xx, domain, fe_values, n_q_points, dofs_per_cell, 0.25, Yleft, Yright, N);
					cout << "IntegralOfTau_xx" << "(0.25) = " << IntegralOfTau_xx << "\n";

					IntegralOfTau_xx = IntegralOfCoordinateX(tau_xx, domain, fe_values, n_q_points, dofs_per_cell, 0.4, Yleft, Yright, N);
					cout << "IntegralOfTau_xx" << "(0.4) = " << IntegralOfTau_xx << "\n";

					IntegralOfTau_xx = IntegralOfCoordinateX(tau_xx, domain, fe_values, n_q_points, dofs_per_cell, 0.5, Yleft, Yright, N);
					cout << "IntegralOfTau_xx" << "(0.5) = " << IntegralOfTau_xx << "\n";




					double IntegralOfTau_xy = 0.0;
					IntegralOfTau_xy = IntegralOfCoordinateX(tau_xy, domain, fe_values, n_q_points, dofs_per_cell, 0.0, Yleft, Yright, N);
					cout << "IntegralOfTau_xy" << "(0.0) = " << IntegralOfTau_xy << "\n";

					IntegralOfTau_xy = IntegralOfCoordinateX(tau_xy, domain, fe_values, n_q_points, dofs_per_cell, 0.166, Yleft, Yright, N);
					cout << "IntegralOfTau_xy" << "(0.166) = " << IntegralOfTau_xy << "\n";

					IntegralOfTau_xy = IntegralOfCoordinateX(tau_xy, domain, fe_values, n_q_points, dofs_per_cell, 0.25, Yleft, Yright, N);
					cout << "IntegralOfTau_xy" << "(0.25) = " << IntegralOfTau_xy << "\n";

					IntegralOfTau_xy = IntegralOfCoordinateX(tau_xy, domain, fe_values, n_q_points, dofs_per_cell, 0.4, Yleft, Yright, N);
					cout << "IntegralOfTau_xy" << "(0.4) = " << IntegralOfTau_xy << "\n";

					IntegralOfTau_xy = IntegralOfCoordinateX(tau_xy, domain, fe_values, n_q_points, dofs_per_cell, 0.5, Yleft, Yright, N);
					cout << "IntegralOfTau_xy" << "(0.5) = " << IntegralOfTau_xy << "\n";



					IntegralOfTau_xy = 0.0;
					IntegralOfTau_xy = IntegralOfCoordinateY(tau_xy, domain, fe_values, n_q_points, dofs_per_cell, 0.0, Xdown, Xup, N);
					cout << "IntegralOfTau_xy(" << 0.0 << ") = " << IntegralOfTau_xy << "\n";

					IntegralOfTau_xy = IntegralOfCoordinateY(tau_xy, domain, fe_values, n_q_points, dofs_per_cell, Yright/4, Xdown, Xup, N);
					cout << "IntegralOfTau_xy(" << Yright/4 << ") = " << IntegralOfTau_xy << "\n";

					IntegralOfTau_xy = IntegralOfCoordinateY(tau_xy, domain, fe_values, n_q_points, dofs_per_cell, Yright/2, Xdown, Xup, N);
					cout << "IntegralOfTau_xy(" << Yright/2 << ") = " << IntegralOfTau_xy << "\n";

					IntegralOfTau_xy = IntegralOfCoordinateY(tau_xy, domain, fe_values, n_q_points, dofs_per_cell, 3*Yright/4, Xdown, Xup, N);
					cout << "IntegralOfTau_xy(" << 3*Yright/4 << ") = " << IntegralOfTau_xy << "\n";

					IntegralOfTau_xy = IntegralOfCoordinateY(tau_xy, domain, fe_values, n_q_points, dofs_per_cell, Yright-0.009, Xdown, Xup, N);
					cout << "IntegralOfTau_xy(" << Yright-0.009 << ") = " << IntegralOfTau_xy << "\n";




					double IntegralOfTau_yy = 0.0;
					IntegralOfTau_yy = IntegralOfCoordinateY(tau_yy, domain, fe_values, n_q_points, dofs_per_cell, 0.0, Xdown, Xup, N);
					cout << "IntegralOfTau_yy(" << 0.0 << ") = " << IntegralOfTau_yy << "\n";

					IntegralOfTau_yy = IntegralOfCoordinateY(tau_yy, domain, fe_values, n_q_points, dofs_per_cell, Yright/4, Xdown, Xup, N);
					cout << "IntegralOfTau_yy(" << Yright/4 << ") = " << IntegralOfTau_yy << "\n";

					IntegralOfTau_yy = IntegralOfCoordinateY(tau_yy, domain, fe_values, n_q_points, dofs_per_cell, Yright/2, Xdown, Xup, N);
					cout << "IntegralOfTau_yy(" << Yright/2 << ") = " << IntegralOfTau_yy << "\n";

					IntegralOfTau_yy = IntegralOfCoordinateY(tau_yy, domain, fe_values, n_q_points, dofs_per_cell, 3*Yright/4, Xdown, Xup, N);
					cout << "IntegralOfTau_yy(" << 3*Yright/4 << ") = " << IntegralOfTau_yy << "\n";

					IntegralOfTau_yy = IntegralOfCoordinateY(tau_yy, domain, fe_values, n_q_points, dofs_per_cell, Yright-0.009, Xdown, Xup, N);
					cout << "IntegralOfTau_yy(" << Yright-0.009 << ") = " << IntegralOfTau_yy << "\n";

				}





		    };

		


		};

		
	return EXIT_SUCCESS;
	}
}
