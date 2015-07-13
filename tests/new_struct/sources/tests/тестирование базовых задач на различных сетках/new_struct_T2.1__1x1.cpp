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
/*
extern void make_grid(
        dealii::Triangulation< 2 >&,
        vec<prmt::Point<2>>,
        vec<st>);

extern void set_grid(
        dealii::Triangulation< 2 >&,
        vec<prmt::Point<2>>,
        vec<prmt::Point<2>>);
*/
void debputs()
{
    static int n = 0;
    printf("DEBUG %d\n", n);
    n++;
};


// template <typename Func, typename Func2>
st foo(cst i, lmbd<st(cst)> &&func, lmbd<st(cst)> &&func2)
{
    return func2(func(i));
};


dbl Analytic_function (const dealii::Point<2> p, cst n, dbl nu)	// Коли
{
    cdbl PI = 3.14159265359;
    dbl Uz = 0.0;
    dbl Uw = 0.0;
    dbl b = 3.0;
    dbl c0 = 0.5;
    dbl C0 = 1.0/3.0;
//    dbl nu = 0.25;
    for (st i = 1; i < n+1; ++i)
    {
        Uw += (nu * b * 4.0 / (std::pow(PI, 3.0) * std::pow((2.0 * i - 1.0), 3.0)) *
                cosh((2.0 * i - 1.0) * PI * p(1)) / sinh((2.0 * i - 1.0) * PI / 2.0 * b) +
        8.0 / (std::pow(PI, 4.0) * std::pow((2.0 * i - 1.0), 4.0))) *
        cos((2.0 * i - 1.0) * PI * p(0));
    };
    Uz = Uw - nu * (  (-std::pow(p(1), 2.0) / 2.0 + C0 ) * (p(0) - c0)  +  std::pow(p(0) - c0, 3.0) / 6.0  );
    // printf("Uber %ld %f %f\n", n, Uw, Uz);
//	Uz = -p(0);
    return Uz;
};










int main()
{
    enum {x, y, z};
    debputs();
    lmbd<st(cst)> add_i = [](cst i){return i-1;};
    printf("%ld\n", foo(10, [](cst i){return i-1;}, [](cst i){return i+3;}));


    //NIKOLA_ELASSTIC_PROBLEM
    if (1)
    {
        Domain<2> domain;
		double Yung		= 1.00;		//-Yung * p(0)
		double Puasson	= 0.25;
        if (1)
		{
			dealii::GridIn<2> gridin;
			gridin.attach_triangulation(domain.grid);
//			std::ifstream f("1x5.msh");
//			std::ifstream f("1x10.msh");
//			std::ifstream f("1x1.msh");
//			std::ifstream f("1x1_T1.2.msh");
//			std::ifstream f("1x2_T1.2.msh");
//			std::ifstream f("1x3_T1.2.msh");
//			std::ifstream f("1x1_T1.2_nestr.msh");
//			std::ifstream f("1x2_T1.2_nestr.msh");
			std::ifstream f("1x1_T2.1.msh");
//			std::ifstream f("test_T2.2.msh");
//			std::ifstream f("10x10.msh");
//			std::ifstream f("2x2.msh");
			Assert (dim==2, ExcInternalError());
			gridin.read_msh(f);
			domain.grid.refine_global(4);
//			domain.grid.refine_global(1);
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
        element_matrix.C .resize (3);
        EPTools ::set_isotropic_elascity{yung : Yung, puasson : Puasson}(element_matrix.C[0]);
//        EPTools ::set_isotropic_elascity{yung : 100.0, puasson : 0.1}(element_matrix.C[1]);

        // T2.2
        vec<arr<typename Nikola::SourceVector<2>::Func, 2>> U(1);
        U[0][x] = [&element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[0][x][x][z][z];}; //Uz
        U[0][y] = [&element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[0][y][y][z][z];}; //Uz
//        U[1][x] = [&element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[1][x][x][z][z];}; //Uz
//        U[1][y] = [&element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[1][y][y][z][z];}; //Uz
        vec<arr<typename Nikola::SourceVector<2>::Func, 2>> tau(1);
        tau[0][x] = [] (const dealii::Point<2> &p) {return 0.0;};	//tau_zx
        tau[0][y] = [] (const dealii::Point<2> &p) {return 0.0;};	//tau_zy
//        tau[1][y] = [] (const dealii::Point<2> &p) {return 0.0;};	//tau_zy
//        tau[1][x] = [] (const dealii::Point<2> &p) {return 0.0;};	//tau_zx

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

		FILE *Ftau_xx, *Ftau_yy, *Ftau_xy, *FU_x_grad, *FU_y_grad, *FU_x, *FU_y, *Ftau_zz;
		FU_x = fopen ("out/FU_x.gpl", "w");					//U_x
		FU_y = fopen ("out/FU_y.gpl", "w");					//U_y
		FU_x_grad = fopen ("out/FU_x_grad.gpl", "w");		//U_x_grad
		FU_y_grad = fopen ("out/FU_y_grad.gpl", "w");		//U_y_grad
		Ftau_xx = fopen ("out/Ftau_xx.gpl", "w");			//tau_xx
		Ftau_yy = fopen ("out/Ftau_yy.gpl", "w");			//tau_yy
		Ftau_xy = fopen ("out/Ftau_xy.gpl", "w");			//tau_xy
		Ftau_zz = fopen ("out/Ftau_zz.gpl", "w");			//tau_zz

				FILE *FAtau_xx, *FAtau_yy, *FAtau_xy, *FAU_x_grad, *FAU_y_grad, *FAU_x, *FAU_y, *FAtau_zz;
				FAU_x = fopen ("out_analytic/FAU_x.gpl", "w");					//U_x
				FAU_y = fopen ("out_analytic/FAU_y.gpl", "w");					//U_y
				FAU_x_grad = fopen ("out_analytic/FAU_x_grad.gpl", "w");		//U_x_grad
				FAU_y_grad = fopen ("out_analytic/FAU_y_grad.gpl", "w");		//U_y_grad
				FAtau_xx = fopen ("out_analytic/FAtau_xx.gpl", "w");			//tau_xx
				FAtau_yy = fopen ("out_analytic/FAtau_yy.gpl", "w");			//tau_yy
				FAtau_xy = fopen ("out_analytic/FAtau_xy.gpl", "w");			//tau_xy
				FAtau_zz = fopen ("out_analytic/FAtau_zz.gpl", "w");			//tau_zz


		int quadrature_points = 0;												//Общее количество квадратурных точек
        dbl integral_x = 0.0;
        dbl integral_y = 0.0;
		dbl Square = 0.0;
		int ID = -100;															//id вершины
		int n_Square = 0;

        dealii::Vector<dbl> x_coordinate_of_vertex(slae.solution.size() / 2);		//массив, содержит координаты x для каждой вершины сетки * 2 (т.к. 2 формулы U_x и U_y)
        dealii::Vector<dbl> y_coordinate_of_vertex(slae.solution.size() / 2);		//массив, содержит координаты y для каждой вершины сетки * 2 (т.к. 2 формулы U_x и U_y)
        dealii::Vector<dbl> tau_xx(slae.solution.size() / 2);						//массив, содержит значения tau_xx
        dealii::Vector<dbl> tau_yy(slae.solution.size() / 2);						//массив, содержит значения tau_yy
        dealii::Vector<dbl> tau_xy(slae.solution.size() / 2);						//массив, содержит значения tau_xy
        dealii::Vector<dbl> tau_zz(slae.solution.size() / 2);						//массив, содержит значения tau_xy
//        dealii::Vector<dbl> test_length(slae.solution.size());					//удалить //для тестирования правильного вывода в файл

//        dealii::Vector<dbl> grad_x(slae.solution.size());						//хранит значения производных по x
//		grad_x = 0.0;
//       dealii::Vector<dbl> grad_y(slae.solution.size());						//хранит значения производных по y
//		grad_y = 0.0;
        dealii::Vector<dbl> U_x_gradX(slae.solution.size() / 2);				//хранит значения производных по x
        dealii::Vector<dbl> U_y_gradX(slae.solution.size() / 2);				//хранит значения производных по x
        dealii::Vector<dbl> U_x_gradY(slae.solution.size() / 2);				//хранит значения производных по x
        dealii::Vector<dbl> U_y_gradY(slae.solution.size() / 2);				//хранит значения производных по x

        dealii::Vector<dbl> U_x(slae.solution.size() / 2);						//хранит значения U_x
        dealii::Vector<dbl> U_y(slae.solution.size() / 2);						//хранит значения U_y


//        dealii::Vector<dbl> s_values_x(slae.solution.size());
//        dealii::Vector<dbl> s_values_y(slae.solution.size());
//        s_values_x = 0.0;
//        s_values_y = 0.0;

        {
            dealii::QGauss<2>  quadrature_formula(2);

            dealii::FEValues<2> fe_values (domain.dof_handler.get_fe(), quadrature_formula,
                    dealii::update_quadrature_points | dealii::update_gradients | dealii::update_JxW_values |
                    dealii::update_values);

            cst dofs_per_cell = fe.dofs_per_cell;								// 8
            cst n_q_points = quadrature_formula.size();							// 4, возвращает кол-во квадратурных точек
//	        dealii::Vector<dbl> cell_grad_x (dofs_per_cell);					//хранит сумму всех значений производной на ячейке
//	        dealii::Vector<dbl> cell_grad_y (dofs_per_cell);					//хранит сумму всех значений производной на ячейке
//	        dealii::Vector<dbl> x_coordinate_of_quadrature_point(n_q_points*domain.grid.n_active_cells());	//массив, содержит координаты x для каждой quad-вершины сетки * 2 (т.к. 2 формулы U_x и U_y)
//	        dealii::Vector<dbl> y_coordinate_of_quadrature_point(n_q_points*domain.grid.n_active_cells());	//массив, содержит координаты y для каждой quad-вершины сетки * 2 (т.к. 2 формулы U_x и U_y)


//			grad_y = 0.0;

            vec<dealii::types::global_dof_index> local_dof_indices (dofs_per_cell);

			printf("\t\tdofs_per_cell = %d\n", dofs_per_cell);
			printf("\t\tn_q_points = %d\n", n_q_points);
			printf("\t\tslae.solution.size() = %d\n", slae.solution.size());

/*
            for (; cell != endc; ++cell)
            {
                fe_values .reinit (cell);
				cell_grad_x = 0.0;
				cell_grad_y = 0.0;
				for (st i = 0; i < dofs_per_cell; ++i)							//	dofs_per_cell = 8,	
					for (st q_point = 0; q_point < n_q_points; ++q_point)		//	n_q_points = 4,	квадратурная формула
						{
							cell_grad_x[q_point] += fe_values.shape_grad (i, q_point)[0] * fe_values.JxW(q_point);
							cell_grad_y[q_point] += fe_values.shape_grad (i, q_point)[1] * fe_values.JxW(q_point);
						}
                cell->get_dof_indices (local_dof_indices);
				for (st i = 0; i < dofs_per_cell; ++i)
					{
	                    grad_x(local_dof_indices[i]) += cell_grad_x(i);
	                    grad_y(local_dof_indices[i]) += cell_grad_y(i);
					}
            };
*/

			int number_of_dots = 0;
            auto cell = domain.dof_handler.begin_active();
            auto endc = domain.dof_handler.end();
			//Пробежать по всем ячейкам, узнать номер DoF, 
			cell = domain.dof_handler.begin_active();
    	    endc = domain.dof_handler.end();
			printf("\t\tnumber_of_dots = %d\n", number_of_dots);




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
			printf("\t\tnumber_of_dots = %d\n", number_of_dots);


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
				fprintf(  TestFile1, "%d\t%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\r\n",	i, x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), slae.solution(i), cos(x_coordinate_of_vertex(i)), -sin(y_coordinate_of_vertex(i)) );
//				fprintf(  TestFile1, "%d\t%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\r\n",	i, x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), slae.solution(i), 2*x_coordinate_of_vertex(i), 2*y_coordinate_of_vertex(i) );
//				fprintf(  TestFile1, "%d\t%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\r\n",	i, x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), slae.solution(i), 				1.0, 					1.0 );
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
			printf("\t\tnumber_of_dots = %d\n", number_of_dots);

//------------------------------------------------------------------------------

/*
			int DD = 0;
			DD = 0 / 2;
			std::cout << 0 << "  0 / 2 = " << DD << "\n";
			DD = 1 / 2;
			std::cout << 0 << "  1 / 2 = " << DD << "\n";
			DD = 2 / 2;
			std::cout << 1 << "  2 / 2 = " << DD << "\n";
			DD = 3 / 2;
			std::cout << 1 << "  3 / 2 = " << DD << "\n";
			DD = 4 / 2;
			std::cout << 2 << "  4 / 2 = " << DD << "\n";
			DD = 5 / 2;
			std::cout << 2 << "  5 / 2 = " << DD << "\n";
			DD = 6 / 2;
			std::cout << 3 << "  6 / 2 = " << DD << "\n";
			DD = 7 / 2;
			std::cout << 3 << "  7 / 2 = " << DD << "\n";
*/



/*
//цикл вывода по всем квадратурным точкам
			for(int i=0; i<n_q_points*domain.grid.n_active_cells(); ++i)
				{
					fprintf(  FgradX,"%9.4f\t%9.4f\t%10.6f\r\n", x_coordinate_of_quadrature_point(i), y_coordinate_of_quadrature_point(i), grad_x(i) );
					fprintf(  FgradY,"%9.4f\t%9.4f\t%10.6f\r\n", x_coordinate_of_quadrature_point(i), y_coordinate_of_quadrature_point(i), grad_y(i) );
				}
*/







/*
			FILE *SolutionVector;												//Вывод вектора решений как есть, с чередованием U_x и U_y
			SolutionVector = fopen ("out/SolutionVector.gpl", "w");
			for(int i=0; i<domain.dof_handler.n_dofs(); ++i)
			{
//																		"x"						"y"					 		"F"			  "NOV"
//				fprintf(  SolutionVector,"%9.4f\t%9.4f\t%9.4f\t%9d\r\n",	x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), slae.solution(i),		i	);
				!(i % 2) ? U_x(i / 2) = slae.solution(i)						//чётное
						 : U_y(i / 2) = slae.solution(i);						//нечётное
			}
			fclose(SolutionVector);


			FILE *outFile1, *outFile2;											//Вывод по отдельности функций U_x и U_y
			outFile1 = fopen ("out/U_x.gpl", "w");
			outFile2 = fopen ("out/U_y.gpl", "w");
			for(int i=0; i<domain.dof_handler.n_dofs() / 2; ++i)				//for(int i=0; i<slae.solution.size() / 2; ++i)
			{
//																		"x"						"y"							"F" 	 "NOV"
				fprintf(  outFile1, "%9.4f\t%9.4f\t%9.4f\t%4d\r\n",	x_coordinate_of_vertex(i), y_coordinate_of_vertex(i),	U_x(i),		i	);		//чётное
				fprintf(  outFile2, "%9.4f\t%9.4f\t%9.4f\t%4d\r\n",	x_coordinate_of_vertex(i), y_coordinate_of_vertex(i),	U_y(i),		i	);		//нечётное
//				cout << "\nU_y = " << U_y(i);
			}
			fclose(outFile1);
			fclose(outFile2);
*/







//------------------------------------------------------------------------------
//Создание файла GRAD.gpl и 
//создать файл GRAD.gpl - содержит производные
//			HCPTools::print_heat_gradient<2>(slae.solution, domain, "GRAD");
//			std::cout << "HCPTools::print_heat_gradient<2> \n";
			HCPTools::print_elasticity_gradient<2>(slae.solution, domain, "GRAD");			//Вызов градиент-функции для U_x
			std::cout << "HCPTools::print_elasticity_gradient<2> \n";
//			GRADIENT::gradient(slae.solution, domain);
//			std::cout << "GRADIENT::gradient\n";


//ввод значений производных из файла GRAD.gpl
//ввод новых значений:
//							x_coordinate_of_vertex(slae.solution.size() / 2)
//							y_coordinate_of_vertex(slae.solution.size() / 2)
//							U_x(slae.solution.size() / 2)
//							U_y(slae.solution.size() / 2)
//							U_x_gradX(slae.solution.size() / 2)
//							U_x_gradY(slae.solution.size() / 2)
//							U_y_gradX(slae.solution.size() / 2)
//							U_y_gradY(slae.solution.size() / 2)
			int w = 0;																	//w - для результата работы функции read_files()
			int M = slae.solution.size(), N = 6;										//M - количество строк в файле, N - количество столбцов в файле. M = 2178, N = 5;
//			w = READ_FILE_TO_ARRAY::read_files(grad_x, grad_y, N);						//передача через указатель. СОБСТВЕННО САМО ЧТЕНИЕ
			w = READ_FILE_TO_ARRAY::read_files(	x_coordinate_of_vertex, y_coordinate_of_vertex,
												U_x, U_y,
												U_x_gradX, U_x_gradY,
												U_y_gradX, U_y_gradY, N);				//передача через указатель. СОБСТВЕННО САМО ЧТЕНИЕ
//------------------------------------------------------------------------------




//------------------------------------------------------------------------------
//Проверка правильности записи в массивы из файла
/*
			std::cout << "\n";
			for(int i=0; i<slae.solution.size() / 2; ++i)
				{
						std::cout << "X (" << i << ") = " << x_coordinate_of_vertex(i) << "\t";
						std::cout << "Y (" << i << ") = " << y_coordinate_of_vertex(i) << "\t";
						std::cout << "U_x  (" << i << ") = " << U_x(i) << "\t\t";
						std::cout << "U_y  (" << i << ") = " << U_y(i) << "\t\t";
						std::cout << "U_x_gradX (" << i << ") = " << U_x_gradX(i) << "\t\t";
						std::cout << "U_y_gradX (" << i << ") = " << U_y_gradX(i) << "\t\t";
						std::cout << "U_x_gradY (" << i << ") = " << U_x_gradY(i) << "\t\t";
						std::cout << "U_y_gradY (" << i << ") = " << U_y_gradY(i) << "\t\t";
						std::cout << "\n";
				}
*/
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

				fprintf( FU_x, "%d%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), U_x(i) );
				fprintf( FU_y, "%d%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), U_y(i) );

//-------------------------------------

				fprintf( FU_x_grad, "%d\t%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), U_x(i), U_x_gradX(i), U_x_gradY(i) );
				fprintf( FU_y_grad, "%d\t%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), U_y(i), U_y_gradX(i), U_y_gradY(i) );

//-------------------------------------

				tau_xx(i)	=	element_matrix.C[0][x][x][x][x] * U_x_gradX(i)
							+	element_matrix.C[0][x][x][y][y] * U_y_gradY(i)			//Здесь grad_y(i) так как длина grad_y равна количеству вершин (slae.solution.size())
							+	element_matrix.C[0][x][x][z][z] * U[0][x](point);

				tau_yy(i)	=	element_matrix.C[0][y][y][x][x] * U_x_gradX(i)			//Здесь grad_y(i) так как длина grad_y равна количеству вершин (slae.solution.size())
							+	element_matrix.C[0][y][y][y][y] * U_y_gradY(i)
							+	element_matrix.C[0][y][y][z][z] * U[0][y](point);
				
				fprintf( Ftau_xx, "%d\t%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), tau_xx(i) );
				fprintf( Ftau_yy, "%d\t%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), tau_yy(i) );

//-------------------------------------

				tau_xy(i)	=	element_matrix.C[0][x][y][x][y] * 
								( U_x_gradY(i) + U_y_gradX(i) );
				fprintf( Ftau_xy, "%d\t%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), tau_xy(i) );

			}
//==============================================================================
//==============================================================================






//==============================================================================
//==============================================================================
//цикл вывода аналитических функций по всем вершинам
			for(int i=0; i<slae.solution.size() / 2; ++i)
			{
//Для того чтобы подставлять значения кординат точек
				dealii::Point<2, double> point(2.0,444.4);
				point(0) = x_coordinate_of_vertex(i);
				point(1) = y_coordinate_of_vertex(i);

//-------------------------------------

				fprintf( FAU_x, "%d%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), Yung / 2 * ( -point(1)*point(1) + point(0)*point(0) + (1*1*1-1)/12 )    );
				fprintf( FAU_y, "%d%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), Yung * point(0) * point(1) );

//-------------------------------------

				fprintf( FAU_x_grad, "%d\t%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), Yung / 2 * ( -point(1)*point(1) + point(0)*point(0) + (1*1*1-1)/12 ), 
															Yung/2*2*point(0),				Yung/2*(-2)*point(1) );
				fprintf( FAU_y_grad, "%d\t%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), Yung * point(0) * point(1), 
															Yung * point(1),				Yung * point(0) );

//-------------------------------------

				fprintf( FAtau_xx, "%d\t%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), 0.0 );
				fprintf( FAtau_yy, "%d\t%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), 0.0 );
				fprintf( FAtau_xy, "%d\t%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), 0.0 );

//-------------------------------------

				fprintf( FAtau_zz, "%d\t%9.4f\t%9.4f\t%9.4f\r\n", i, point(0), point(1), Puasson );

			}
//==============================================================================
//==============================================================================







/*
			//тот самых цикл вывода по всем вершинам
			for(int i=0; i<slae.solution.size(); ++i)
			{
//Ftest_length -  для тестирования правильного вывода в файл
//					fprintf(  Ftest_length,"%9.4f\t%9.4f\t%9.4f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), test_length(i) );

//					!(i % 2) ? fprintf(  Ftau_xx,"%9.4f\t%9.4f\t%9.4f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), test_length(i) )
//							 : fprintf(  Ftau_yy,"%9.4f\t%9.4f\t%9.4f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), test_length(i) );


//Для того чтобы подставлять значения кординат точек
					dealii::Point<2, double> point(2.0,444.4);
					point(0) = x_coordinate_of_vertex(i);
					point(1) = y_coordinate_of_vertex(i);
//slae.solution(i) - это массив, с чередованием функций {U_x, U_y, U_x, U_y, ..}
//Поэтому далее, в некоторых случаях: 		U_x = slae.solution(i),
//									 		U_y = slae.solution(i+1)

//-------------------------------------

				!(i % 2) ? fprintf(  FgradX,"%9.4f\t%9.4f\t%9.4f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), grad_x(i) )
						 : fprintf(  FgradY,"%9.4f\t%9.4f\t%9.4f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), grad_y(i) );

//-------------------------------------

				!(i % 2) ?
				tau_xx(i)	=	element_matrix.C[0][x][x][x][x] * slae.solution(i)		*	grad_x(i)
							+	element_matrix.C[0][x][x][y][y] * slae.solution(i+1)	*	grad_y(i)			//Здесь grad_y(i) так как длина grad_y равна количеству вершин (slae.solution.size())
							+	element_matrix.C[0][x][x][z][z] * U[0][x](point)
				: tau_xx(i)	=	0.0;

				!(i % 2) ?
				tau_yy(i)	=	0.0 :
				tau_yy(i)	=	element_matrix.C[0][y][y][x][x] * slae.solution(i-1)	*	grad_x(i)			//Здесь grad_y(i) так как длина grad_y равна количеству вершин (slae.solution.size())
							+	element_matrix.C[0][y][y][y][y] * slae.solution(i)		*	grad_y(i)
							+	element_matrix.C[0][y][y][z][z] * U[0][x](point);
				
				!(i % 2) ? fprintf(  Ftau_xx,"%9.4f\t%9.4f\t%9.4f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), tau_xx(i) )
						 : fprintf(  Ftau_yy,"%9.4f\t%9.4f\t%9.4f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), tau_yy(i) );

//-------------------------------------

				!(i % 2) ?
				tau_xy(i)	=	element_matrix.C[0][x][y][x][y] * 
								(
									slae.solution(i)	*	grad_y(i) +
									slae.solution(i+1)	*	grad_x(i)
								)
				: tau_xy(i)	=	0.0;

				!(i % 2) ? fprintf(  Ftau_xy,"%9.4f\t%9.4f\t%10.6f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), tau_xy(i) )
						 : tau_xy(i)	=	0.0;

//-------------------------------------

				!(i % 2) ? fprintf(  FU_x,"%9.4f\t%9.4f\t%9.4f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), slae.solution(i) )
						 : fprintf(  FU_y,"%9.4f\t%9.4f\t%9.4f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), slae.solution(i) );

//-------------------------------------

				!(i % 2) ?
				tau_zz(i)	=	element_matrix.C[0][x][x][z][z] * tau_xx(i) +
								element_matrix.C[0][y][y][z][z] * tau_yy(i) +
								element_matrix.C[0][z][z][z][z] * U[0][x](point)
//					tau_zz(i)	=	nu_xz * tau_xx + nu_yz * tau_yy + E_z * U_z
								
				: tau_zz(i)	=	0.0;

				!(i % 2) ? fprintf(  Ftau_zz,"%9.4f\t%9.4f\t%10.6f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), tau_zz(i) )
						 : tau_zz(i)	=	0.0;
			}
*/


        };

		
//        EPTools::print_move<2> (slae.solution, domain.dof_handler, "move.gpd");
		fclose(FU_x);
		fclose(FU_y);
		fclose(FU_x_grad);
		fclose(FU_y_grad);
		fclose(Ftau_xx);
		fclose(Ftau_yy);
		fclose(Ftau_xy);
		fclose(Ftau_zz);

			fclose(FAU_x);
			fclose(FAU_y);
			fclose(FAU_x_grad);
			fclose(FAU_y_grad);
			fclose(FAtau_xx);
			fclose(FAtau_yy);
			fclose(FAtau_xy);
			fclose(FAtau_zz);


		printf("\t\tquadrature_points = %d\n", quadrature_points);
		printf("\t\tSquare = %f\n", Square);
		printf("\t\tn_Square = %d\n", n_Square);
		printf("\t\tintegral_x = %f\n", integral_x);
		printf("\t\tintegral_y = %f\n", integral_y);

    };

    
return EXIT_SUCCESS;
}
