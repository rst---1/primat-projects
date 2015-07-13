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

namespace NEW_STRUCT_TEMPERATURE_T1_1
{
	double Analytic_functionU_z (const dealii::Point<2> p, cst n, double Puasson, double c0, double C0, double bb)
	{
		cdbl PI = 3.14159265359;
		dbl Uz = 0.0;
		dbl Uw = 0.0;
		for (st i = 1; i < n+1; ++i)
		{
		    Uw += (Puasson * bb * 4.0 / (std::pow(PI, 3.0) * std::pow((2.0 * i - 1.0), 3.0)) *
		            cosh((2.0 * i - 1.0) * PI * p(1)) / sinh((2.0 * i - 1.0) * PI / 2.0 * bb) +
		    8.0 / (std::pow(PI, 4.0) * std::pow((2.0 * i - 1.0), 4.0))) *
		    cos( (2.0 * i - 1.0) * PI * p(0) );
		};
		Uz = Uw - Puasson * (  (-std::pow(p(1), 2.0) / 2.0 + C0 ) * (p(0) - c0)  +  std::pow(p(0) - c0, 3.0) / 6.0  );
		return Uz;
	};

	double Analytic_functionU_z_gradX (const dealii::Point<2> p, cst n, double Puasson, double c0, double C0, double bb)
	{
		cdbl PI = 3.14159265359;
		dbl Uz = 0.0;
		dbl Uw = 0.0;
		for (st i = 1; i < n+1; ++i)
		{
		    Uw += (Puasson * bb * 4.0 / (std::pow(PI, 3.0) * std::pow((2.0 * i - 1.0), 3.0)) *
		            cosh((2.0 * i - 1.0) * PI * p(1)) / sinh((2.0 * i - 1.0) * PI / 2.0 * bb) +
		    8.0 / (std::pow(PI, 4.0) * std::pow((2.0 * i - 1.0), 4.0))) *
		    sin( (2.0 * i - 1.0) * PI * p(0) ) * (-1) * ((2.0 * i - 1.0) * PI);
		};
		Uz = Uw - Puasson * (  (-std::pow(p(1), 2.0) / 2.0 + C0 )  +  3.0 * std::pow(p(0) - c0, 2.0) / 6.0  );
		return Uz;
	};

	double Analytic_functionU_z_gradY (const dealii::Point<2> p, cst n, double Puasson, double c0, double C0, double bb)
	{
		cdbl PI = 3.14159265359;
		dbl Uz = 0.0;
		dbl Uw = 0.0;
		for (st i = 1; i < n+1; ++i)
		{
		    Uw += ( Puasson * bb * 4.0 / (std::pow(PI, 3.0) * std::pow((2.0 * i - 1.0), 3.0)) *
		            ((2.0 * i - 1.0) * PI) * sinh((2.0 * i - 1.0) * PI * p(1)) / sinh((2.0 * i - 1.0) * PI / 2.0 * bb) ) *
		    cos( (2.0 * i - 1.0) * PI * p(0) );
		};
		Uz = Uw - Puasson * (  (-p(1) * 2.0 / 2.0 ) * (p(0) - c0)  );
		return Uz;
	};

	double Analytic_functionTau_zx (const dealii::Point<2> p, cst n, double Puasson, double c0, double C0, double bb, double mu)
	{
		cdbl PI = 3.14159265359;
		dbl Uw = 0.0;
		for (st i = 1; i < n+1; ++i)
		{
		    Uw += (Puasson * bb * 4.0 / (std::pow(PI, 3.0) * std::pow((2.0 * i - 1.0), 3.0)) *
		            cosh((2.0 * i - 1.0) * PI * p(1)) / sinh((2.0 * i - 1.0) * PI / 2.0 * bb) +
		    8.0 / (std::pow(PI, 4.0) * std::pow((2.0 * i - 1.0), 4.0))) *
		    sin( (2.0 * i - 1.0) * PI * p(0) ) * (-1) * ( (2.0 * i - 1.0) * PI );
		};
		return mu*Uw;
	};

	double Analytic_functionTau_zy (const dealii::Point<2> p, cst n, double Puasson, double c0, double C0, double bb, double mu)
	{
		cdbl PI = 3.14159265359;
		dbl Uw = 0.0;
		for (st i = 1; i < n+1; ++i)
		{
		    Uw += ( Puasson * bb * 4.0 / (std::pow(PI, 3.0) * std::pow((2.0 * i - 1.0), 3.0)) *
		            ((2.0 * i - 1.0) * PI) * sinh((2.0 * i - 1.0) * PI * p(1)) / sinh((2.0 * i - 1.0) * PI / 2.0 * bb) ) *
		    cos( (2.0 * i - 1.0) * PI * p(0) );
		};
		Uw += 2 * Puasson * p(1) * (p(0) - c0);
		return mu*Uw;
	};



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
		for (auto cell = domain.dof_handler.begin_active(); cell != domain.dof_handler.end(); ++cell)
		{
			fe_values .reinit (cell);
			cell->get_dof_indices(local_dof_indices);
			for(int i = 0; i < n_q_points; ++i)
			{
				integral_U_z += Vec(local_dof_indices[i])*fe_values.JxW(i);
				Square += fe_values.JxW(i);
			}
		}
		return integral_U_z;
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
		std::cout << "\n\nВызов функции main() из файла new_struct_temperature.cpp" << "\n\n";
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


		//HEAT_CONDUCTION_NIKOLA_PROBLEM
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
				dealii::GridIn<2> gridin;
				gridin.attach_triangulation(domain.grid);
				std::ifstream f(FileOfGrid);
				std::cout << "\tFileOfGrid = " << FileOfGrid << "\n";			//output
				Assert (dim==2, ExcInternalError());
				gridin.read_msh(f);
				domain.grid.refine_global(RefineNum);
					std::cout << "\tNumber of active cells:       "				//output
				        	<< domain.grid.n_active_cells()
				        	<< std::endl;
		    }

		    dealii::FE_Q<2> fe(1);
		    domain.dof_init (fe);

		    SystemsLinearAlgebraicEquations slae;
		    ATools ::trivial_prepare_system_equations (slae, domain);


		    LaplacianScalar<2> element_matrix (domain.dof_handler.get_fe());
		    {
				element_matrix.C .resize(3);										//???????
				element_matrix.C[0][x][x] = mu1;
				element_matrix.C[0][x][y] = 0.0;
				element_matrix.C[0][y][x] = 0.0;
				element_matrix.C[0][y][y] = mu1;
					element_matrix.C[1][x][x] = mu2;
					element_matrix.C[1][x][y] = 0.0;
					element_matrix.C[1][y][x] = 0.0;
					element_matrix.C[1][y][y] = mu2;
				element_matrix.C[2][x][x] = mu3;
				element_matrix.C[2][x][y] = 0.0;
				element_matrix.C[2][y][x] = 0.0;
				element_matrix.C[2][y][y] = mu3;
		    };

			dealii::Vector<dbl> U_x(slae.solution.size());						//Вводим новый вектор
			dealii::Vector<dbl> U_y(slae.solution.size());						//Вводим новый вектор
			dealii::Vector<dbl> tau_zz(slae.solution.size());						//Вводим новый вектор
//			READ_FILE_TO_ARRAY::read_vector("out/FU_x_81.gpl", U_x, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести
//			READ_FILE_TO_ARRAY::read_vector("out/FU_y_81.gpl", U_y, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести
//			READ_FILE_TO_ARRAY::read_vector("out/Ftau_zz_81.gpl", tau_zz, 3, 4);		//Считываем в новый вектор данные из файла. ...3, 6 - четвёртый столбец из шести

		    vec<arr<typename Nikola::SourceScalar<2>::Func, 2>> U(3);
		    U[0][x] = [mu1,  &element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[0][x][x]*1.0;};		//Ux
		    U[0][y] = [mu1,  &element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[0][x][x]*0.0;};		//Uy
		    U[1][x] = [mu2,  &element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[1][x][x]*1.0;};		//Ux
		    U[1][y] = [mu2,  &element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[1][x][x]*0.0;};		//Uy
		    U[2][x] = [mu3,  &element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[2][x][x]*1.0;};		//Ux
		    U[2][y] = [mu3,  &element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[2][x][x]*0.0;};		//Uy
//		    U[0][x] = [&U_x, &domain, &element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[0][x][x] * HCPTools::print_function_temperature<2>(U_x, domain, p);}; //Ux
//		    U[0][y] = [&U_y, &domain, &element_matrix] (const dealii::Point<2> &p) {return element_matrix.C[0][y][y] * HCPTools::print_function_temperature<2>(U_y, domain, p);}; //Uy


		    vec<typename Nikola::SourceScalar<2>::Func> tau(3);
        	tau[0] = [Yung1, c0] (const dealii::Point<2> &p) {return 0.0;};
        	tau[1] = [Yung2, c0] (const dealii::Point<2> &p) {return 0.0;};
        	tau[2] = [Yung3, c0] (const dealii::Point<2> &p) {return 0.0;};


		    Nikola::SourceScalar<2> element_rhsv (U, tau, domain.dof_handler.get_fe());
		    Assembler::assemble_matrix<2> (slae.matrix, element_matrix, domain.dof_handler);
		    Assembler::assemble_rhsv<2> (slae.rhsv, element_rhsv, domain.dof_handler);


		    dealii::SolverControl solver_control (10000, eps);
		    dealii::SolverCG<> solver (solver_control);
		    solver.solve (
		            slae.matrix,
		            slae.solution,
		            slae.rhsv,
		            dealii::PreconditionIdentity()
		            );
			cout << "\teps = " << eps << "\n";									//output
			cout << "\n\n";

	//==============================================================================



			int quadrature_points = 0;												//Общее количество квадратурных точек
			int ID = -100;															//id вершины

		    dealii::Vector<dbl> x_coordinate_of_vertex(slae.solution.size());		//массив, содержит координаты x
		    dealii::Vector<dbl> y_coordinate_of_vertex(slae.solution.size());		//массив, содержит координаты y

	//		dealii::Vector<dbl> U_x(slae.solution.size());						//хранит значения U_x
	//		dealii::Vector<dbl> U_y(slae.solution.size());						//хранит значения U_y
			dealii::Vector<dbl> U_z(slae.solution.size());						//хранит значения U_z
	//        dealii::Vector<dbl> U_x_gradX(slae.solution.size());				//хранит значения производных по x
	//        dealii::Vector<dbl> U_x_gradY(slae.solution.size());				//хранит значения производных по y
	//        dealii::Vector<dbl> U_y_gradX(slae.solution.size());				//хранит значения производных по x
	//        dealii::Vector<dbl> U_y_gradY(slae.solution.size());				//хранит значения производных по y
		    dealii::Vector<dbl> U_z_gradX(slae.solution.size());				//хранит значения производных по x
		    dealii::Vector<dbl> U_z_gradY(slae.solution.size());				//хранит значения производных по y
	//        dealii::Vector<dbl> tau_xx(slae.solution.size());					//массив, содержит значения tau_xx
	//        dealii::Vector<dbl> tau_yy(slae.solution.size());					//массив, содержит значения tau_yy
	//        dealii::Vector<dbl> tau_xy(slae.solution.size());					//массив, содержит значения tau_xy
	//        dealii::Vector<dbl> tau_zz(slae.solution.size());					//массив, содержит значения tau_xy
		    dealii::Vector<dbl> tau_zx(slae.solution.size());					//массив, содержит значения tau_zx
		    dealii::Vector<dbl> tau_zy(slae.solution.size());					//массив, содержит значения tau_zy

					//		dealii::Vector<dbl> AU_x(slae.solution.size());						//хранит значения AU_x
					//		dealii::Vector<dbl> AU_y(slae.solution.size());						//хранит значения AU_y
							dealii::Vector<dbl> AU_z(slae.solution.size());						//хранит значения AU_z
					//        dealii::Vector<dbl> AU_x_gradX(slae.solution.size());				//хранит значения производных по x
					//        dealii::Vector<dbl> AU_x_gradY(slae.solution.size());				//хранит значения производных по y
					//        dealii::Vector<dbl> AU_y_gradX(slae.solution.size());				//хранит значения производных по x
					//        dealii::Vector<dbl> AU_y_gradY(slae.solution.size());				//хранит значения производных по y
							dealii::Vector<dbl> AU_z_gradX(slae.solution.size());				//хранит значения производных по x
							dealii::Vector<dbl> AU_z_gradY(slae.solution.size());				//хранит значения производных по y
					//        dealii::Vector<dbl> Atau_xx(slae.solution.size());				//массив, содержит значения Atau_xx
					//        dealii::Vector<dbl> Atau_yy(slae.solution.size());				//массив, содержит значения Atau_yy
					//        dealii::Vector<dbl> Atau_xy(slae.solution.size());				//массив, содержит значения Atau_xy
					//        dealii::Vector<dbl> Atau_zz(slae.solution.size());				//массив, содержит значения Atau_xy
							dealii::Vector<dbl> Atau_zx(slae.solution.size());					//массив, содержит значения Atau_zx
							dealii::Vector<dbl> Atau_zy(slae.solution.size());					//массив, содержит значения Atau_zy



		    {
		        dealii::QGauss<2>  quadrature_formula(2);

		        dealii::FEValues<2> fe_values (domain.dof_handler.get_fe(), quadrature_formula,
		                dealii::update_quadrature_points | dealii::update_gradients | dealii::update_JxW_values |
		                dealii::update_values);

		        cst dofs_per_cell = fe.dofs_per_cell;								// 8
		        cst n_q_points = quadrature_formula.size();							// 4, возвращает кол-во квадратурных точек

				std::cout << "dofs_per_cell = " << dofs_per_cell << "\n";									//4						//output
				std::cout << "n_q_points = " << n_q_points << "\n";											//4						//output
				std::cout << "domain.dof_handler.n_dofs() = " << domain.dof_handler.n_dofs() << "\n";		//117					//output
				std::cout << "slae.solution.size() = " << slae.solution.size() << "\n";						//117					//output


		        vec<dealii::types::global_dof_index> local_dof_indices (dofs_per_cell);

		        auto cell = domain.dof_handler.begin_active();
		        auto endc = domain.dof_handler.end();
				//Пробежать по всем ячейкам, узнать номер DoF, 
				cell = domain.dof_handler.begin_active();
			    endc = domain.dof_handler.end();





	//------------------------------------------------------------------------------
	//инедексация векторов
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
	//				test_length(ID) = ID;
					x_coordinate_of_vertex(ID) = cell->vertex(0)[0];
					y_coordinate_of_vertex(ID) = cell->vertex(0)[1];
						ID = cell->vertex_dof_index (1, 0);
	//					test_length(ID) = ID;
						x_coordinate_of_vertex(ID) = cell->vertex(1)[0];
						y_coordinate_of_vertex(ID) = cell->vertex(1)[1];
							ID = cell->vertex_dof_index (2, 0);
	//						test_length(ID) = ID;
							x_coordinate_of_vertex(ID) = cell->vertex(2)[0];
							y_coordinate_of_vertex(ID) = cell->vertex(2)[1];
								ID = cell->vertex_dof_index (3, 0);
	//							test_length(ID) = ID;
								x_coordinate_of_vertex(ID) = cell->vertex(3)[0];
								y_coordinate_of_vertex(ID) = cell->vertex(3)[1];

				}

	//------------------------------------------------------------------------------











	//------------------------------------------------------------------------------
	//Создание файла GRAD.gpl и 
	//создать файл GRAD.gpl - содержит производные
				HCPTools::print_heat_gradient<2>(slae.solution, domain, "GRAD");
				std::cout << "HCPTools::print_heat_gradient<2> \n";

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
				w = READ_FILE_TO_ARRAY::read_files_for_temperature(	x_coordinate_of_vertex, y_coordinate_of_vertex,
																U_z, U_z_gradX, U_z_gradY, N);				//передача через указатель. СОБСТВЕННО САМО ЧТЕНИЕ
	//------------------------------------------------------------------------------










	//==============================================================================
	//==============================================================================
	//цикл вывода аналитических функций по всем вершинам
				for(int i=0; i<slae.solution.size(); ++i)
				{
	//Для того чтобы подставлять значения кординат точек
					dealii::Point<2, double> point(2.0,444.4);
					point(0) = x_coordinate_of_vertex(i);
					point(1) = y_coordinate_of_vertex(i);

	//-------------------------------------

					AU_z(i) = -(point(0)-c0);
					fprintf( FAU_z, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), AU_z(i)  );

	//-------------------------------------

					AU_z_gradX(i) = -1.0;
					AU_z_gradY(i) = 0.0;
					fprintf( FAU_z_grad, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1),
															AU_z(i), AU_z_gradX(i), AU_z_gradY(i) );

	//-------------------------------------

					Atau_zx(i) = 0.0;
					Atau_zy(i) = 0.0;
					fprintf( FAtau_zx, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), Atau_zx(i) );
					fprintf( FAtau_zy, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), Atau_zy(i) );

				}

	//==============================================================================
	//==============================================================================







	//------------------------------------------------------------------------------
	//Вычисление интеграла
			{
				double AddConstant = 0.0;
				double integralConst = 0.0;
				cout << " => (in beginning) Integral(AU_z) = " << Integral(AU_z, domain, fe_values, n_q_points, dofs_per_cell) << "\n";

				for(int j = 0; j < 15; j++)
				{
					integralConst = 0.95 * Integral(AU_z, domain, fe_values, n_q_points, dofs_per_cell) / bb;
					Add_to_vector(AU_z, -integralConst);
					AddConstant += -integralConst;
				}
				integralConst = Integral(AU_z, domain, fe_values, n_q_points, dofs_per_cell);
				cout << " => Integral(AU_z) = " << integralConst << "\n";
				cout << "AddConstant = " << AddConstant << "\n";
				std::cout << "\n";
			}
	//------------------------------------------------------------------------------


/*
	//------------------------------------------------------------------------------
	//Вычисление интеграла
				{
					std::cout << "\n";
					double	integralConst = Integral(U_z, domain, fe_values, n_q_points, dofs_per_cell);
					int ID = 0;
					int iter = 0;
					int iterOFwhile = 0;
					double AveArea = (Xup - Xdown)*(Yright - Yleft) / U_z.size(); 				//Найдём среднюю площадь на одну вершину. Для этого общуую площадь разделим на кол-во вершин.
					double rx = sqrt(AveArea) / 50;
					double ry = sqrt(AveArea) / 50;
					double Xcentr = (Xup + Xdown) / 2;											//Найдём центр 
					double Ycentr = (Yright + Yleft) / 2;										//Найдём центр 
					dealii::Point<2, double> point(2.0,444.4);

							cout << "AveArea = " << AveArea << "\n";																		//output
							cout << "rx = " << rx << "\n";																					//output
							cout << "ry = " << ry << "\n";																					//output
							ID = ID_reserch(U_z, x_coordinate_of_vertex, y_coordinate_of_vertex, Xcentr, Ycentr, rx, ry, point);
							while(ID == 0)
							{
								iter++;
								rx *= 1.05; ry *= 1.05;
								ID = ID_reserch(U_z, x_coordinate_of_vertex, y_coordinate_of_vertex, Xcentr, Ycentr, rx, ry, point);
							}

							cout << " x_coordinate_of_vertex(ID) = " <<  x_coordinate_of_vertex(ID) << "\n";								//output
							cout << " y_coordinate_of_vertex(ID) = " <<  y_coordinate_of_vertex(ID) << "\n";								//output
							cout << " Xcentr = " <<  Xcentr << "\n";																		//output
							cout << " Ycentr = " <<  Ycentr << "\n";																		//output
							cout << "iter = " << iter << "\n";																				//output
							cout << "ID = " << ID << "\n";																					//output
							cout << "U_z(ID) = " << U_z(ID) << "\n";
							cout << "AU_z(ID) = " << AU_z(ID) << "\n";
							cout << "U_z(ID) - AU_z(ID) = " << U_z(ID) - AU_z(ID) << "\n";
							std::cout << "\n";

								integralConst = 0;
								while(  ( abs(integralConst) > 1000.0 ) && (iterOFwhile < 50)  )
								{
									iterOFwhile++;
									if (integralConst > 0.0) { Add_to_vector(U_z, -100.0); }
										else { Add_to_vector(U_z, 100.0); }
									integralConst = Integral(U_z, domain, fe_values, n_q_points, dofs_per_cell);																	//output
									cout << " => Integral = " << integralConst << "\n";
									std::cout << "\n";
								}

											integralConst = U_z(ID) - AU_z(ID);
											cout << "integralConst = U_z(ID) - AU_z(ID) = " << integralConst << "\n";
											cout << "ID = " << ID << "\n";																					//output
											cout << "U_z(ID) = " << U_z(ID) << "\n";	
											cout << "AU_z(ID) = " << AU_z(ID) << "\n";
											Add_to_vector(U_z, -integralConst);
											integralConst = Integral(U_z, domain, fe_values, n_q_points, dofs_per_cell);																	//output
											cout << " => Integral(U_z) = " << integralConst << "\n";
				//							integralConst = U_z(ID) - AU_z(ID);
											std::cout << "\n";
																//output
				}
	//------------------------------------------------------------------------------
*/

	//------------------------------------------------------------------------------
	//Вычисление интеграла
			{
				double AddConstant = 0.0;
				double integralConst = 0.0;
				cout << " => (in beginning) Integral(U_z) = " << Integral(U_z, domain, fe_values, n_q_points, dofs_per_cell) << "\n";

				for(int j = 0; j < 15; j++)
				{
					integralConst = 0.95 * Integral(U_z, domain, fe_values, n_q_points, dofs_per_cell) / bb;
					Add_to_vector(U_z, -integralConst);
					AddConstant += -integralConst;
				}
				integralConst = Integral(U_z, domain, fe_values, n_q_points, dofs_per_cell);	
				cout << " => Integral(U_z) = " << integralConst << "\n";
				cout << "AddConstant = " << AddConstant << "\n";
				std::cout << "\n";
			}
	//------------------------------------------------------------------------------









	//==============================================================================
	//==============================================================================
	double		deltABSU_z = -100;
	double		deltAVEU_z = 0.0;
	double		deltOTNU_z = -100;
	double		deltMaxABSU_z = -100;
	double		deltMaxOTNU_z = -100;

	double		deltABSU_z_gradX = -100;
	double		deltOTNU_z_gradX = -100;
	double		deltMaxABSU_z_gradX = -100;
	double		deltMaxOTNU_z_gradX = -100;

	double		deltABSU_z_gradY = -100;
	double		deltOTNU_z_gradY = -100;
	double		deltMaxABSU_z_gradY = -100;
	double		deltMaxOTNU_z_gradY = -100;

	double		deltABStau_zx = -100;
	double		deltAVEtau_zx = 0.0;
	double		deltOTNtau_zx = -100;
	double		deltMaxABStau_zx = -100;
	double		deltMaxOTNtau_zx = -100;

	double		deltABStau_zy = -100;
	double		deltAVEtau_zy = 0.0;
	double		deltOTNtau_zy = -100;
	double		deltMaxABStau_zy = -100;
	double		deltMaxOTNtau_zy = -100;

	//тот самых цикл вывода по всем вершинам
				for(int i=0; i<slae.solution.size(); ++i)
				{
	//Для того чтобы подставлять значения кординат точек
					dealii::Point<2, double> point(2.0,444.4);
					point(0) = x_coordinate_of_vertex(i);
					point(1) = y_coordinate_of_vertex(i);

	//-------------------------------------

					fprintf( FU_z, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), U_z(i) );		//U_z
					
	//-------------------------------------

					fprintf( FU_z_grad, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), U_z(i), U_z_gradX(i), U_z_gradY(i) );
					if( (point(0)>Xdown)&&(point(0)<Xup) && (point(1)>Yleft)&&(point(1)<Yright) )
					{
						deltABSU_z = abs(AU_z(i)-U_z(i));
						deltAVEU_z += deltABSU_z;
						if(U_z(i) != 0.0) {deltOTNU_z = abs(deltABSU_z / U_z(i) * 100.0);}	else deltOTNU_z = 0.0;
						if (deltMaxABSU_z < deltABSU_z) {deltMaxABSU_z = deltABSU_z;}
						if (deltMaxOTNU_z < deltOTNU_z) {deltMaxOTNU_z = deltOTNU_z;}

						deltABSU_z_gradX = abs(AU_z_gradX(i)-U_z_gradX(i));
						if(U_z_gradX(i) != 0.0) {deltOTNU_z_gradX = abs(deltABSU_z_gradX / U_z_gradX(i) * 100.0);}	else deltOTNU_z_gradX = 0.0;
						if (deltMaxABSU_z_gradX < deltABSU_z_gradX) {deltMaxABSU_z_gradX = deltABSU_z_gradX;}
						if (deltMaxOTNU_z_gradX < deltOTNU_z_gradX) {deltMaxOTNU_z_gradX = deltOTNU_z_gradX;}

						deltABSU_z_gradY = abs(AU_z_gradY(i)-U_z_gradY(i));
						if(U_z_gradY(i) != 0.0) {deltOTNU_z_gradY = abs(deltABSU_z_gradY / U_z_gradY(i) * 100.0);}	else deltOTNU_z_gradY = 0.0;
						if (deltMaxABSU_z_gradY < deltABSU_z_gradY) {deltMaxABSU_z_gradY = deltABSU_z_gradY;}
						if (deltMaxOTNU_z_gradY < deltOTNU_z_gradY) {deltMaxOTNU_z_gradY = deltOTNU_z_gradY;}
					}

	//-------------------------------------

					if(point(0) < h)
					{
						tau_zx(i)	=	element_matrix.C[0][x][x] * ( U[0][x](point) / element_matrix.C[0][x][x] + U_z_gradX(i) );	//Разделим, т. к. при вызове функции возвращается значение,
																																	//умноженное на element_matrix.C[0][x][x]
						tau_zy(i)	=	element_matrix.C[0][x][x] * ( U[0][y](point) / element_matrix.C[0][y][y] + U_z_gradY(i) );	//Разделим, т. к. при вызове функции возвращается значение,
																																	//умноженное на element_matrix.C[0][y][y]
					}
					else if(point(0) < 2*h)
					{
						tau_zx(i)	=	element_matrix.C[1][x][x] * ( U[1][x](point) / element_matrix.C[1][x][x] + U_z_gradX(i) );	//Разделим, т. к. при вызове функции возвращается значение,
																																	//умноженное на element_matrix.C[0][x][x]
						tau_zy(i)	=	element_matrix.C[1][x][x] * ( U[1][y](point) / element_matrix.C[1][y][y] + U_z_gradY(i) );	//Разделим, т. к. при вызове функции возвращается значение,
																																	//умноженное на element_matrix.C[0][y][y]
					}
					else
					{
						tau_zx(i)	=	element_matrix.C[2][x][x] * ( U[2][x](point) / element_matrix.C[2][x][x] + U_z_gradX(i) );	//Разделим, т. к. при вызове функции возвращается значение,
																																	//умноженное на element_matrix.C[0][x][x]
						tau_zy(i)	=	element_matrix.C[2][x][x] * ( U[2][y](point) / element_matrix.C[2][y][y] + U_z_gradY(i) );	//Разделим, т. к. при вызове функции возвращается значение,
																																	//умноженное на element_matrix.C[0][y][y]
					}
				
					fprintf( Ftau_zx, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), tau_zx(i) );
					fprintf( Ftau_zy, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), tau_zy(i) );
					if( (point(0)>Xdown)&&(point(0)<Xup) && (point(1)>Yleft)&&(point(1)<Yright) )
					{
						deltABStau_zx = abs(Atau_zx(i)-tau_zx(i));
						deltAVEtau_zx += deltABStau_zx;
						if(tau_zx(i) != 0.0) {deltOTNtau_zx = abs(deltABStau_zx / tau_zx(i) * 100.0);}	else deltOTNtau_zx = 0.0;
						if (deltMaxABStau_zx < deltABStau_zx) {deltMaxABStau_zx = deltABStau_zx;}
						if (deltMaxOTNtau_zx < deltOTNtau_zx) {deltMaxOTNtau_zx = deltOTNtau_zx;}

						deltABStau_zy = abs(Atau_zy(i)-tau_zy(i));
						deltAVEtau_zy += deltABStau_zy;
						if(tau_zy(i) != 0.0) {deltOTNtau_zy = abs(deltABStau_zy / tau_zy(i) * 100.0);}	else deltOTNtau_zy = 0.0;
						if (deltMaxABStau_zy < deltABStau_zy) {deltMaxABStau_zy = deltABStau_zy;}
						if (deltMaxOTNtau_zy < deltOTNtau_zy) {deltMaxOTNtau_zy = deltOTNtau_zy;}
					}

				}
				std::cout << "deltMaxABSU_z = " << deltMaxABSU_z << "\n";										//output
				std::cout << "deltAVEU_z = " << deltAVEU_z / slae.solution.size() << "\n";						//output
				std::cout << "deltMaxOTNU_z = " << deltMaxOTNU_z << "%\n";										//output

				std::cout << "deltMaxABSU_z_gradX = " << deltMaxABSU_z_gradX << "\n";							//output
				std::cout << "deltMaxOTNU_z_gradX = " << deltMaxOTNU_z_gradX << "%\n";							//output

				std::cout << "deltMaxABSU_z_gradY = " << deltMaxABSU_z_gradY << "\n";							//output
				std::cout << "deltMaxOTNU_z_gradY = " << deltMaxOTNU_z_gradY << "%\n";							//output

				std::cout << "deltMaxABStau_zx = " << deltMaxABStau_zx << "\n";									//output
				std::cout << "deltAVEtau_zx = " << deltAVEtau_zx / slae.solution.size() << "\n";				//output
				std::cout << "deltMaxOTNtau_zx = " << deltMaxOTNtau_zx << "%\n";								//output

				std::cout << "deltMaxABStau_zy = " << deltMaxABStau_zy << "\n";									//output
				std::cout << "deltAVEtau_zy = " << deltAVEtau_zy / slae.solution.size() << "\n";				//output
				std::cout << "deltMaxOTNtau_zy = " << deltMaxOTNtau_zy << "%\n";								//output

				std::cout << "\tFileOfGrid			= " << FileOfGrid << "_" << domain.grid.n_active_cells() << "_temperture" <<"\n";
	//==============================================================================
	//==============================================================================






	//==============================================================================
	//==============================================================================
/*
		FILE *ErrFU_z;				ErrFU_z = fopen ("out_errors/ErrFU_z.gpl", "w+");
		FILE *ErrFtau_zx;			ErrFtau_zx = fopen ("out_errors/ErrFtau_zx.gpl", "w+");
		FILE *ErrFtau_zy;			ErrFtau_zy = fopen ("out_errors/ErrFtau_zy.gpl", "w+");

	double		deltABS = -100;
	double		deltOTN = -100;
	//цикл вывода аналитических функций по всем вершинам
				for(int i=0; i<slae.solution.size(); ++i)
				{
	//Для того чтобы подставлять значения кординат точек
					dealii::Point<2, double> point(2.0,444.4);
					point(0) = x_coordinate_of_vertex(i);
					point(1) = y_coordinate_of_vertex(i);

	//-------------------------------------

					deltABS = abs(AU_z(i)-U_z(i));
					if(AU_z(i) != 0.0) {deltOTN = abs(deltABS / AU_z(i) * 100.0);}
					else deltOTN = 0.0;
					fprintf( ErrFU_z, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), U_z(i), AU_z(i), deltABS, deltOTN  );

	//-------------------------------------

					deltABS = abs(Atau_zx(i)-tau_zx(i));
					if(Atau_zx(i) != 0.0) {deltOTN = abs(deltABS / Atau_zx(i) * 100.0);}
					else deltOTN = 0.0;
					fprintf( ErrFtau_zx, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), tau_zx(i), Atau_zx(i), deltABS, deltOTN );

					deltABS = abs(Atau_zy(i)-tau_zy(i));
					if(Atau_zy(i) != 0.0) {deltOTN = abs(deltABS / Atau_zy(i) * 100.0);}
					else deltOTN = 0.0;
					fprintf( ErrFtau_zy, "%8d\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\t\t%5.35f\r\n", i, point(0), point(1), tau_zy(i), Atau_zy(i), deltABS, deltOTN );

				}
			fclose(ErrFU_z);
			fclose(ErrFtau_zx);
			fclose(ErrFtau_zy);
*/
	//==============================================================================
	//==============================================================================


		    };

		


	//		printf("\t\tquadrature_points = %d\n", quadrature_points);
	//		printf("\t\tSquare = %f\n", Square);
	//		printf("\t\tn_Square = %d\n", n_Square);
	//		printf("\t\tintegral_x = %f\n", integral_x);
	//		printf("\t\tintegral_y = %f\n", integral_y);

		};

    
	return EXIT_SUCCESS;
	}
}
