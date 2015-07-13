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


namespace NEW_STRUCT_ELASTIC
{
	void Add_to_vector(dealii::Vector<dbl> &Vec, const double number)						//Позволяет добавить к решению константу
	{
		for (int i=0; i<Vec.size(); ++i)	
		{
			Vec(i) += number;
		}
	}

	void Integral(dealii::Vector<dbl> &Vec, Domain<2> &domain, dealii::FEValues<2> &fe_values, cst &n_q_points, cst &dofs_per_cell)						//Позволяет добавить к решению константу
	{
		vec<dealii::types::global_dof_index> local_dof_indices (dofs_per_cell);
		double integral_U_z = 0.0;
		for (auto cell = domain.dof_handler.begin_active(); cell != domain.dof_handler.end(); ++cell)
		{
			fe_values .reinit (cell);
			cell->get_dof_indices(local_dof_indices);
			for(int i = 0; i < n_q_points; ++i)
			{
				integral_U_z += Vec(local_dof_indices[i])*fe_values.JxW(i);
			}
		}
		std::cout << "\n";
		std::cout << "U_z = " << integral_U_z << "\n";
		std::cout << "\n";
	}





	int main(		const double &Yung,
					const double &Puasson,
					const double &c0,
					const double &C0,
					const double &bb,
					const double &mu,
					string FileOfGrid,
					const double &Yleft,
					const double &Yright,
					const double &Xup,
					const double &Xdown,
					int &RefineNum,
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
				dealii::GridIn<2> gridin;
				gridin.attach_triangulation(domain.grid);							//INPUT
				std::ifstream f(FileOfGrid);
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
		    element_matrix.C .resize (1);											//один материал
		    EPTools ::set_isotropic_elascity{yung : Yung, puasson : Puasson}(element_matrix.C[0]);
	//        EPTools ::set_isotropic_elascity{yung : 100.0, puasson : 0.1}(element_matrix.C[1]);


			cout << "\n\n";
//			dealii::Vector<dbl> U_z(slae.solution.size() / 2);						//Вводим новый вектор
			dealii::Vector<dbl> U_z(slae.solution.size() / 2);						//Вводим новый вектор
			dealii::Vector<dbl> tau_zx(slae.solution.size() / 2);						//Вводим новый вектор
			dealii::Vector<dbl> tau_zy(slae.solution.size() / 2);						//Вводим новый вектор
//			READ_FILE_TO_ARRAY::read_vector("out/FU_z_4225.gpl", U_z, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести
//			READ_FILE_TO_ARRAY::read_vector("out/Ftau_zx_4225.gpl", tau_zx, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести
//			READ_FILE_TO_ARRAY::read_vector("out/Ftau_zy_4225.gpl", tau_zy, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести

			READ_FILE_TO_ARRAY::read_vector("out/FU_z.gpl", U_z, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести
			READ_FILE_TO_ARRAY::read_vector("out/Ftau_zx.gpl", tau_zx, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести
			READ_FILE_TO_ARRAY::read_vector("out/Ftau_zy.gpl", tau_zy, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести

		    vec<arr<typename Nikola::SourceVector<2>::Func, 2>> U(1);
//		    U[0][x] = [&U_z, &domain, &element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[0][x][x][z][z] * HCPTools::print_function_elastic<2>(U_z, domain, p);}; //Uz					//INPUT
//		    U[0][y] = [&U_z, &domain, &element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[0][x][x][z][z] * HCPTools::print_function_elastic<2>(U_z, domain, p);}; //Uz					//INPUT
		    U[0][x] = [&element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[0][x][x][z][z] * (-1)*(p(0)-c0);}; //Uz					//INPUT
		    U[0][y] = [&element_matrix, c0] (const dealii::Point<2> &p) {return element_matrix.C[0][y][y][z][z] * (-1)*(p(0)-c0);}; //Uz					//INPUT

		    vec<arr<typename Nikola::SourceVector<2>::Func, 2>> tau(1);
//		    tau[0][x] = [&tau_zx, &domain] (const dealii::Point<2> &p) {return HCPTools::print_function_elastic<2>(tau_zx, domain, p);};				//tau_zx
//		    tau[0][y] = [&tau_zy, &domain] (const dealii::Point<2> &p) {return HCPTools::print_function_elastic<2>(tau_zy, domain, p);};				//tau_zy
		    tau[0][x] = [] (const dealii::Point<2> &p) {return 0.0;};	//tau_zx													//INPUT
		    tau[0][y] = [] (const dealii::Point<2> &p) {return 0.0;};	//tau_zy													//INPUT

		    Nikola::SourceVector<2> element_rhsv (U, tau, domain.dof_handler.get_fe());

		    Assembler ::assemble_matrix<2> (slae.matrix, element_matrix, domain.dof_handler);
		    Assembler ::assemble_rhsv<2> (slae.rhsv, element_rhsv, domain.dof_handler);

		    dealii::SolverControl solver_control (10000, 1e-12);
		    dealii::SolverCG<> solver (solver_control);
		    solver.solve (
		            slae.matrix,
		            slae.solution,
		            slae.rhsv
		            ,dealii::PreconditionIdentity()
		            );
	//==============================================================================



			int quadrature_points = 0;												//Общее количество квадратурных точек
//		    dbl integral_x = 0.0;
//		    dbl integral_y = 0.0;
//			dbl Square = 0.0;
			int ID = -100;															//id вершины
//			int n_Square = 0;

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
//				printf("\t\tnumber_of_dots = %8d\n", number_of_dots);



	//------------------------------------------------------------------------------
	//Вычисление интеграла
				{
					Integral(slae.solution, domain, fe_values, n_q_points, dofs_per_cell);

					Add_to_vector(slae.solution, IntegralCoefficient);

					Integral(slae.solution, domain, fe_values, n_q_points, dofs_per_cell);
				}
	//------------------------------------------------------------------------------


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
				std::cout << "READ_FILE_TO_ARRAY::read_files_for_elelastic \n";
					//Прибавление константы к функции, чтобы приравнять интеграл функции к нулю
	//				for(int i=0; i<slae.solution.size() / 2; ++i)							//INPUT
	//				{
	//					U_x(i) -= 0.125;
	//				}
	//------------------------------------------------------------------------------
			    











	//==============================================================================
	//==============================================================================
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
					if ( (i == 10) ) {cout << "U_x() = " << U_x(i) << "\n";}
	//-------------------------------------

					fprintf( FU_x_grad, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), U_x(i), U_x_gradX(i), U_x_gradY(i) );
					fprintf( FU_y_grad, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), U_y(i), U_y_gradX(i), U_y_gradY(i) );

	//-------------------------------------			//INPUT

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
//								+	element_matrix.C[0][y][y][z][z] * (-1) * (point(0) - c0);
				
					fprintf( Ftau_xx, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), tau_xx(i) );
					fprintf( Ftau_yy, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), tau_yy(i) );

	//-------------------------------------			//INPUT

					tau_xy(i)	=	element_matrix.C[0][x][y][x][y] * 
									( U_x_gradY(i) + U_y_gradX(i) );
					fprintf( Ftau_xy, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), tau_xy(i) );

	//-------------------------------------			//INPUT

					tau_zz(i)	=	element_matrix.C[0][x][z][x][z] * tau_xx(i)
								+	element_matrix.C[0][y][z][y][z] * tau_yy(i)
//								+	Yung * (-1) * (point(0) - c0);
								+	Yung * U[0][x](point) / element_matrix.C[0][x][x][z][z];	//Разделим, т. к. при вызове функции возвращается значение,
																								//умноженное на element_matrix.C[0][x][x][z][z]
	//							+	element_matrix.C[0][z][z][z][z] * (-1) * (point(0) - c0);
					fprintf( Ftau_zz, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1),  tau_zz(i));

				}
				std::cout << "Output files \n";
	//==============================================================================
	//==============================================================================














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



	//цикл вывода аналитических функций по всем вершинам
				for(int i=0; i<slae.solution.size() / 2; ++i)
				{
	//Для того чтобы подставлять значения кординат точек
					dealii::Point<2, double> point(2.0,444.4);
					point(0) = x_coordinate_of_vertex(i);
					point(1) = y_coordinate_of_vertex(i);

	//-------------------------------------			//INPUT

					AU_x(i) = Puasson / 2 * ( -point(1)*point(1) + (point(0) - c0)*(point(0) - c0) + (bb*bb*bb - 1.0)/12 );
					fprintf( FAU_x, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), AU_x(i) );


					AU_y(i) = Puasson * (point(0) - c0) * point(1);
					fprintf( FAU_y, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), AU_y(i) );

	//-------------------------------------			//INPUT

					AU_x_gradX(i) = Puasson / 2 * 2 * (point(0) - c0);
					AU_x_gradY(i) = Puasson / 2 * (-2) * point(1);
					fprintf( FAU_x_grad, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1),
																AU_x(i), AU_x_gradX(i), AU_x_gradY(i) );
					if( (point(0)>Xdown)&&(point(0)<Xup) && (point(1)>Yleft)&&(point(1)<Yright) )
					{
						deltABSU_x = abs(AU_x(i)-U_x(i));
						if(U_x(i) != 0.0) {deltOTNU_x = abs(deltABSU_x / U_x(i) * 100.0);}	else deltOTNU_x = 0.0;
						if (deltMaxABSU_x < deltABSU_x) {deltMaxABSU_x = deltABSU_x;}
						if (deltMaxOTNU_x < deltOTNU_x) {deltMaxOTNU_x = deltOTNU_x;}

						deltABSU_x_gradX = abs(AU_x_gradX(i)-U_x_gradX(i));
						if(U_x_gradX(i) != 0.0) {deltOTNU_x_gradX = abs(deltABSU_x_gradX / U_x_gradX(i) * 100.0);}	else deltOTNU_x_gradX = 0.0;
						if (deltMaxABSU_x_gradX < deltABSU_x_gradX) {deltMaxABSU_x_gradX = deltABSU_x_gradX;}
						if (deltMaxOTNU_x_gradX < deltOTNU_x_gradX) {deltMaxOTNU_x_gradX = deltOTNU_x_gradX;}

						deltABSU_x_gradY = abs(AU_x_gradY(i)-U_x_gradY(i));
						if(U_x_gradY(i) != 0.0) {deltOTNU_x_gradY = abs(deltABSU_x_gradY / U_x_gradY(i) * 100.0);}	else deltOTNU_x_gradY = 0.0;
						if (deltMaxABSU_x_gradY < deltABSU_x_gradY) {deltMaxABSU_x_gradY = deltABSU_x_gradY;}
						if (deltMaxOTNU_x_gradY < deltOTNU_x_gradY) {deltMaxOTNU_x_gradY = deltOTNU_x_gradY;}
					}


					AU_y_gradX(i) = Puasson * point(1);
					AU_y_gradY(i) = Puasson * (point(0) - c0);
					fprintf( FAU_y_grad, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1),
																AU_y(i), AU_y_gradX(i), AU_y_gradY(i) );
					if( (point(0)>Xdown)&&(point(0)<Xup) && (point(1)>Yleft)&&(point(1)<Yright) )
					{
						deltABSU_y = abs(AU_y(i)-U_y(i));
						if(U_y(i) != 0.0) {deltOTNU_y = abs(deltABSU_y / U_y(i) * 100.0);}	else deltOTNU_y = 0.0;
						if (deltMaxABSU_y < deltABSU_y) {deltMaxABSU_y = deltABSU_y;}
						if (deltMaxOTNU_y < deltOTNU_y) {deltMaxOTNU_y = deltOTNU_y;}

						deltABSU_y_gradX = abs(AU_y_gradX(i)-U_y_gradX(i));
						if(U_y_gradX(i) != 0.0) {deltOTNU_y_gradX = abs(deltABSU_y_gradX / U_y_gradX(i) * 100.0);}	else deltOTNU_y_gradX = 0.0;
						if (deltMaxABSU_y_gradX < deltABSU_y_gradX) {deltMaxABSU_y_gradX = deltABSU_y_gradX;}
						if (deltMaxOTNU_y_gradX < deltOTNU_y_gradX) {deltMaxOTNU_y_gradX = deltOTNU_y_gradX;}

						deltABSU_y_gradY = abs(AU_y_gradY(i)-U_y_gradY(i));
						if(U_y_gradY(i) != 0.0) {deltOTNU_y_gradY = abs(deltABSU_y_gradY / U_y_gradY(i) * 100.0);}	else deltOTNU_y_gradY = 0.0;
						if (deltMaxABSU_y_gradY < deltABSU_y_gradY) {deltMaxABSU_y_gradY = deltABSU_y_gradY;}
						if (deltMaxOTNU_y_gradY < deltOTNU_y_gradY) {deltMaxOTNU_y_gradY = deltOTNU_y_gradY;}
					}

	//-------------------------------------

					Atau_xx(i) = 0.0;
					fprintf( FAtau_xx, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), Atau_xx(i) );

					Atau_yy(i) = 0.0;
					fprintf( FAtau_yy, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), Atau_yy(i) );

					Atau_xy(i) = 0.0;
					fprintf( FAtau_xy, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), Atau_xy(i) );

					if( (point(0)>Xdown)&&(point(0)<Xup) && (point(1)>Yleft)&&(point(1)<Yright) )
					{
						deltABStau_xx = abs(Atau_xx(i)-tau_xx(i));
						deltAVEtau_xx += deltABStau_xx;
						if(tau_xx(i) != 0.0) {deltOTNtau_xx = abs(deltABStau_xx / tau_xx(i) * 100.0);}	else deltOTNtau_xx = 0.0;
						if (deltMaxABStau_xx < deltABStau_xx) {deltMaxABStau_xx = deltABStau_xx;}
						if (deltMaxOTNtau_xx < deltOTNtau_xx) {deltMaxOTNtau_xx = deltOTNtau_xx;}

						deltABStau_yy = abs(Atau_yy(i)-tau_yy(i));
						deltAVEtau_yy += deltABStau_yy;
						if(tau_yy(i) != 0.0) {deltOTNtau_yy = abs(deltABStau_yy / tau_yy(i) * 100.0);}	else deltOTNtau_yy = 0.0;
						if (deltMaxABStau_yy < deltABStau_yy) {deltMaxABStau_yy = deltABStau_yy;}
						if (deltMaxOTNtau_yy < deltOTNtau_yy) {deltMaxOTNtau_yy = deltOTNtau_yy;}

						deltABStau_xy = abs(Atau_xy(i)-tau_xy(i));
						deltAVEtau_xy += deltABStau_xy;
						if(tau_xy(i) != 0.0) {deltOTNtau_xy = abs(deltABStau_xy / tau_xy(i) * 100.0);}	else deltOTNtau_xy = 0.0;
						if (deltMaxABStau_xy < deltABStau_xy) {deltMaxABStau_xy = deltABStau_xy;}
						if (deltMaxOTNtau_xy < deltOTNtau_xy) {deltMaxOTNtau_xy = deltOTNtau_xy;}
					}

	//-------------------------------------			//INPUT

					Atau_zz(i) = -Yung*(point(0) - c0);
					fprintf( FAtau_zz, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), Atau_zz(i) );

					if( (point(0)>Xdown)&&(point(0)<Xup) && (point(1)>Yleft)&&(point(1)<Yright) )
					{
						deltABStau_zz = abs(Atau_zz(i)-tau_zz(i));
						deltAVEtau_zz += deltABStau_zz;
						if(tau_zz(i) != 0.0) {deltOTNtau_zz = abs(deltABStau_zz / tau_zz(i) * 100.0);}	else deltOTNtau_zz = 0.0;
						if (deltMaxABStau_zz < deltABStau_zz) {deltMaxABStau_zz = deltABStau_zz;}
						if (deltMaxOTNtau_zz < deltOTNtau_zz) {deltMaxOTNtau_zz = deltOTNtau_zz;}
					}
				}
				std::cout << "\n";
				std::cout << "\n";
				std::cout << "deltMaxABSU_x = " << deltMaxABSU_x << "\n";
				std::cout << "deltMaxOTNU_x = " << deltMaxOTNU_x << "%\n";
				std::cout << "\n";

				std::cout << "deltMaxABSU_x_gradX = " << deltMaxABSU_x_gradX << "\n";
				std::cout << "deltMaxOTNU_x_gradX = " << deltMaxOTNU_x_gradX << "%\n";
				std::cout << "\n";

				std::cout << "deltMaxABSU_x_gradY = " << deltMaxABSU_x_gradY << "\n";
				std::cout << "deltMaxOTNU_x_gradY = " << deltMaxOTNU_x_gradY << "%\n";
				std::cout << "\n";

				std::cout << "deltMaxABSU_y = " << deltMaxABSU_y << "\n";
				std::cout << "deltMaxOTNU_y = " << deltMaxOTNU_y << "%\n";
				std::cout << "\n";

				std::cout << "deltMaxABSU_y_gradX = " << deltMaxABSU_y_gradX << "\n";
				std::cout << "deltMaxOTNU_y_gradX = " << deltMaxOTNU_y_gradX << "%\n";
				std::cout << "\n";

				std::cout << "deltMaxABSU_y_gradY = " << deltMaxABSU_y_gradY << "\n";
				std::cout << "deltMaxOTNU_y_gradY = " << deltMaxOTNU_y_gradY << "%\n";
				std::cout << "\n";

				std::cout << "deltMaxABStau_xx = " << deltMaxABStau_xx << "\n";
				std::cout << "deltAVEtau_xx = " << deltAVEtau_xx / (slae.solution.size() / 2) << "\n";
				std::cout << "deltMaxOTNtau_xx = " << deltMaxOTNtau_xx << "%\n";
				std::cout << "\n";

				std::cout << "deltMaxABStau_yy = " << deltMaxABStau_yy << "\n";
				std::cout << "deltAVEtau_yy = " << deltAVEtau_yy / (slae.solution.size() / 2) << "\n";
				std::cout << "deltMaxOTNtau_yy = " << deltMaxOTNtau_yy << "%\n";
				std::cout << "\n";

				std::cout << "deltMaxABStau_xy = " << deltMaxABStau_xy << "\n";
				std::cout << "deltAVEtau_xy = " << deltAVEtau_xy / (slae.solution.size() / 2) << "\n";
				std::cout << "deltMaxOTNtau_xy = " << deltMaxOTNtau_xy << "%\n";
				std::cout << "\n";

				std::cout << "deltMaxABStau_zz = " << deltMaxABStau_zz << "\n";
				std::cout << "deltAVEtau_zz = " << deltAVEtau_zz / (slae.solution.size() / 2) << "\n";
				std::cout << "deltMaxOTNtau_zz = " << deltMaxOTNtau_zz << "%\n";
				std::cout << "\n";

	//==============================================================================
	//==============================================================================


		    };

		


		};

		
	return EXIT_SUCCESS;
	}
}
