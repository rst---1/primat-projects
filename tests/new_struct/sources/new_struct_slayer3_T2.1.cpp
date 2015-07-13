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
//#include "new_struct0.cpp"
//#include "new_struct_elastic.cpp"
#include "new_struct_elastic_T2.1_slayer3.cpp"									//INPUT
//#include "new_struct_elastic_T2.2_slayer3.cpp"									//INPUT

//#include "new_struct_temperature.cpp"
//#include "new_struct_temperature_T1.1_slayer3.cpp"									//INPUT
#include "new_struct_temperature_T1.2_slayer3.cpp"									//INPUT
//#include "new_struct_temperature_T1.1.cpp"
//#include "new_struct_temperature_T1.2.cpp"


int main1(	int RefineNum, const double bb, const double Yleft, const double Yright,
			const double Xdown, const double Xup, const string FileOfGrid,
			const double Yung1, const double Yung2, const double Yung3,
			const double c0, const double C0, const double eps, const char* FILEBEGIN )
{
	enum {x, y, z};
	std::cout << "\n\n================================================================================\n";
	std::cout << "================================================================================\n";

		const double Puasson1	= 0.4;		//nu									//INPUT
		const double Puasson2	= 0.1;		//nu									//INPUT
		const double Puasson3	= 0.4;		//nu									//INPUT

	    const double mu1 = Yung1 / ( 2.0 * ( 1.0 + Puasson1 ) );
	    const double mu2 = Yung2 / ( 2.0 * ( 1.0 + Puasson2 ) );
	    const double mu3 = Yung3 / ( 2.0 * ( 1.0 + Puasson3 ) );
	    const double h = 0.33;				//толщина слоя							//INPUT

		double IntegralCoefficient = 0.0;										//На сколько приподнять функцию U_z



		char* FILENAME = new char[150];				//содержит начальное имя файла
		strcpy(FILENAME, FILEBEGIN);
		strcat(FILENAME, "FU_x.gpl");
		std::cout << "\tFILENAME = " << FILENAME << "\n";


		FILE *FU_x;			strcpy(FILENAME, FILEBEGIN);		strcat(FILENAME, "FU_x.gpl");							FU_x = fopen (FILENAME, "w+");					//U_x
		FILE *FU_y;			strcpy(FILENAME, FILEBEGIN);		strcat(FILENAME, "FU_y.gpl");							FU_y = fopen (FILENAME, "w+");					//U_y
		FILE *FU_z;			strcpy(FILENAME, FILEBEGIN);		strcat(FILENAME, "FU_z.gpl");							FU_z = fopen (FILENAME, "w+");					//U_z
		FILE *FU_x_grad;			strcpy(FILENAME, FILEBEGIN);		strcat(FILENAME, "FU_x_grad.gpl");							FU_x_grad = fopen (FILENAME, "w+");					//U_x_grad
		FILE *FU_y_grad;			strcpy(FILENAME, FILEBEGIN);		strcat(FILENAME, "FU_y_grad.gpl");							FU_y_grad = fopen (FILENAME, "w+");					//U_y_grad
		FILE *FU_z_grad;			strcpy(FILENAME, FILEBEGIN);		strcat(FILENAME, "FU_z_grad.gpl");							FU_z_grad = fopen (FILENAME, "w+");					//U_z_grad
		FILE *Ftau_xx;			strcpy(FILENAME, FILEBEGIN);		strcat(FILENAME, "Ftau_xx.gpl");							Ftau_xx = fopen (FILENAME, "w+");					//tau_xx
		FILE *Ftau_yy;			strcpy(FILENAME, FILEBEGIN);		strcat(FILENAME, "Ftau_yy.gpl");							Ftau_yy = fopen (FILENAME, "w+");					//tau_yy
		FILE *Ftau_xy;			strcpy(FILENAME, FILEBEGIN);		strcat(FILENAME, "Ftau_xy.gpl");							Ftau_xy = fopen (FILENAME, "w+");					//tau_xy
		FILE *Ftau_zz;			strcpy(FILENAME, FILEBEGIN);		strcat(FILENAME, "Ftau_zz.gpl");							Ftau_zz = fopen (FILENAME, "w+");					//tau_zz
		FILE *Ftau_zx;			strcpy(FILENAME, FILEBEGIN);		strcat(FILENAME, "Ftau_zx.gpl");							Ftau_zx = fopen (FILENAME, "w+");					//tau_zx
		FILE *Ftau_zy;			strcpy(FILENAME, FILEBEGIN);		strcat(FILENAME, "Ftau_zy.gpl");							Ftau_zy = fopen (FILENAME, "w+");					//tau_zy
/*
		FILE *FU_x;				FU_x = fopen ("out/FU_x.gpl", "w+");					//U_y
		FILE *FU_y;				FU_y = fopen ("out/FU_y.gpl", "w+");					//U_y
		FILE *FU_z;				FU_z = fopen ("out/FU_z.gpl", "w+");					//U_z
		FILE *FU_x_grad;		FU_x_grad = fopen ("out/FU_x_grad.gpl", "w+");		//U_x_grad
		FILE *FU_y_grad;		FU_y_grad = fopen ("out/FU_y_grad.gpl", "w+");		//U_y_grad
		FILE *FU_z_grad;		FU_z_grad = fopen ("out/FU_z_grad.gpl", "w+");		//U_z_grad
		FILE *Ftau_xx;			Ftau_xx = fopen ("out/Ftau_xx.gpl", "w+");			//tau_xx
		FILE *Ftau_yy;			Ftau_yy = fopen ("out/Ftau_yy.gpl", "w+");			//tau_yy
		FILE *Ftau_xy;			Ftau_xy = fopen ("out/Ftau_xy.gpl", "w+");			//tau_xy
		FILE *Ftau_zz;			Ftau_zz = fopen ("out/Ftau_zz.gpl", "w+");			//tau_zz
		FILE *Ftau_zx;			Ftau_zx = fopen ("out/Ftau_zx.gpl", "w+");			//tau_zx
		FILE *Ftau_zy;			Ftau_zy = fopen ("out/Ftau_zy.gpl", "w+");			//tau_zy
*/
				FILE *FAU_x;			FAU_x = fopen ("out_analytic/FAU_x.gpl", "w+");					//U_x
				FILE *FAU_y;			FAU_y = fopen ("out_analytic/FAU_y.gpl", "w+");					//U_y
				FILE *FAU_z;			FAU_z = fopen ("out_analytic/FAU_z.gpl", "w+");					//U_z
				FILE *FAU_x_grad;		FAU_x_grad = fopen ("out_analytic/FAU_x_grad.gpl", "w+");		//U_x_grad
				FILE *FAU_y_grad;		FAU_y_grad = fopen ("out_analytic/FAU_y_grad.gpl", "w+");		//U_y_grad
				FILE *FAU_z_grad;		FAU_z_grad = fopen ("out_analytic/FAU_z_grad.gpl", "w+");		//U_z_grad
				FILE *FAtau_xx;			FAtau_xx = fopen ("out_analytic/FAtau_xx.gpl", "w+");			//tau_xx
				FILE *FAtau_yy;			FAtau_yy = fopen ("out_analytic/FAtau_yy.gpl", "w+");			//tau_yy
				FILE *FAtau_xy;			FAtau_xy = fopen ("out_analytic/FAtau_xy.gpl", "w+");			//tau_xy
				FILE *FAtau_zz;			FAtau_zz = fopen ("out_analytic/FAtau_zz.gpl", "w+");			//tau_zz
				FILE *FAtau_zx;			FAtau_zx = fopen ("out_analytic/FAtau_zx.gpl", "w+");			//tau_zx
				FILE *FAtau_zy;			FAtau_zy = fopen ("out_analytic/FAtau_zy.gpl", "w+");			//tau_zy

		FILE *FILES[] = {FU_x, FU_y, FU_z, FU_x_grad, FU_y_grad, FU_z_grad,
						Ftau_xx, Ftau_yy, Ftau_xy, Ftau_zz, Ftau_zx, Ftau_zy,
						FAU_x, FAU_y, FAU_z, FAU_x_grad, FAU_y_grad, FAU_z_grad,
						FAtau_xx, FAtau_yy, FAtau_xy, FAtau_zz, FAtau_zx, FAtau_zy};



//TEMPERATURE BLOCK
	if (0)									//INPUT
	{
		//Зададим пока что функции U[][] и tau[][] через формулы, а не через файлы
		NEW_STRUCT_TEMPERATURE::main(Yung1, Yung2, Yung3, Puasson1, Puasson2, Puasson3, c0, C0, bb, mu1, mu2, mu3, h, FileOfGrid, Yleft, Yright, Xup, Xdown, RefineNum, IntegralCoefficient, eps, FILES);
	}
//end of TEMPERATURE BLOCK

//Нужно закрыть следующие файлы, чтобы можно было корректно из них читать данные
			fclose(FU_z);
			fclose(Ftau_zx);
			fclose(Ftau_zy);

//ELASTIC BLOCK
	if (1)									//INPUT
	{
		NEW_STRUCT_ELASTIC::main(Yung1, Yung2, Yung3, Puasson1, Puasson2, Puasson3, c0, C0, bb, mu1, mu2, mu3, h, FileOfGrid, Yleft, Yright, Xup, Xdown, RefineNum, IntegralCoefficient, eps, FILES);
	}
//end of ELASTIC BLOCK



			fclose(FU_x);
			fclose(FU_y);
//			fclose(FU_z);
			fclose(FU_x_grad);
			fclose(FU_y_grad);
			fclose(FU_z_grad);
			fclose(Ftau_xx);
			fclose(Ftau_yy);
			fclose(Ftau_xy);
			fclose(Ftau_zz);
//			fclose(Ftau_zx);
//			fclose(Ftau_zy);
				fclose(FAU_x);
				fclose(FAU_y);
				fclose(FAU_z);
				fclose(FAU_x_grad);
				fclose(FAU_y_grad);
				fclose(FAU_z_grad);
				fclose(FAtau_xx);
				fclose(FAtau_yy);
				fclose(FAtau_xy);
				fclose(FAtau_zz);
				fclose(FAtau_zx);
				fclose(FAtau_zy);

	
	std::cout << "\n\n================================================================================\n";
	std::cout << "================================================================================\n";
	return EXIT_SUCCESS;
}

void main2(const double Yung1, const double Yung2, const double Yung3)
{
	double bb;
	double Yleft;
	double Yright;
	double Xdown;
	double Xup;
	string FileOfGrid;
	double c0;
	double C0;
	double eps;
		char *FILEBEGIN = new char[100];
		char *STRBEGIN1;
		char *STRBEGIN2;
		char *STRBEGIN3;







	STRBEGIN1 = "out/";
	STRBEGIN2 = "T2.1/";							//INPUT







	STRBEGIN3 = "1x10_unregular/";
	strcpy(FILEBEGIN, STRBEGIN1);					//Запись пути файла
	strcat(FILEBEGIN, STRBEGIN2);					//Запись пути файла
	strcat(FILEBEGIN, STRBEGIN3);					//Запись пути файла
	bb = 10.0;	Yleft = -5.0; Yright = 5.0; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x10_slayer3_unregular.msh"; c0 = 0.5; C0 = 4.325; eps = 1e-12;							//INPUT
	main1( 1, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps, FILEBEGIN );


	STRBEGIN3 = "1x5_unregular/";
	strcpy(FILEBEGIN, STRBEGIN1);					//Запись пути файла
	strcat(FILEBEGIN, STRBEGIN2);					//Запись пути файла
	strcat(FILEBEGIN, STRBEGIN3);					//Запись пути файла
	bb = 5.0;	Yleft = -2.5; Yright = 2.5; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x5_slayer3_unregular.msh"; c0 = 0.5; C0 = 263.0/240.0; eps = 1e-12;							//INPUT
	main1( 1, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps, FILEBEGIN );


	STRBEGIN3 = "1x1_unregular/";
	strcpy(FILEBEGIN, STRBEGIN1);					//Запись пути файла
	strcat(FILEBEGIN, STRBEGIN2);					//Запись пути файла
	strcat(FILEBEGIN, STRBEGIN3);					//Запись пути файла
	bb = 1.0;	Yleft = -0.5; Yright = 0.5; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x1_slayer3_unregular.msh"; c0 = 0.5; C0 = 0.0; eps = 1e-12;							//INPUT
	main1( 1, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps, FILEBEGIN );


	STRBEGIN3 = "1x10/";
	strcpy(FILEBEGIN, STRBEGIN1);					//Запись пути файла
	strcat(FILEBEGIN, STRBEGIN2);					//Запись пути файла
	strcat(FILEBEGIN, STRBEGIN3);					//Запись пути файла
	bb = 10.0;	Yleft = -5.0; Yright = 5.0; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x10_slayer3.msh"; c0 = 0.5; C0 = 4.325; eps = 1e-12;							//INPUT
	main1( 5, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps, FILEBEGIN );		//14


	STRBEGIN3 = "1x5/";
	strcpy(FILEBEGIN, STRBEGIN1);					//Запись пути файла
	strcat(FILEBEGIN, STRBEGIN2);					//Запись пути файла
	strcat(FILEBEGIN, STRBEGIN3);					//Запись пути файла
	bb = 5.0;	Yleft = -2.5; Yright = 2.5; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x5_slayer3.msh"; c0 = 0.5; C0 = 263.0/240.0; eps = 1e-12;							//INPUT
	main1( 5, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps, FILEBEGIN );		//10


	STRBEGIN3 = "1x1/";
	strcpy(FILEBEGIN, STRBEGIN1);					//Запись пути файла
	strcat(FILEBEGIN, STRBEGIN2);					//Запись пути файла
	strcat(FILEBEGIN, STRBEGIN3);					//Запись пути файла
	bb = 1.0;	Yleft = -0.5; Yright = 0.5; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x1_slayer3.msh"; c0 = 0.5; C0 = 0.0; eps = 1e-12;							//INPUT
	main1( 5, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps, FILEBEGIN );		//5




/*
	bb = 10.0;	Yleft = -5.0; Yright = 5.0; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x10_slayer3_unregular.msh"; c0 = 0.5; C0 = 4.325; eps = 1e-12;							//INPUT
	main1( 1, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );

	bb = 5.0;	Yleft = -2.5; Yright = 2.5; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x5_slayer3_unregular.msh"; c0 = 0.5; C0 = 263.0/240.0; eps = 1e-12;							//INPUT
	main1( 1, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );

	bb = 1.0;	Yleft = -0.5; Yright = 0.5; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x1_slayer3_unregular.msh"; c0 = 0.5; C0 = 0.0; eps = 1e-12;							//INPUT
	main1( 1, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );


	bb = 10.0;	Yleft = -5.0; Yright = 5.0; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x10_slayer3.msh"; c0 = 0.5; C0 = 4.325; eps = 1e-12;							//INPUT
	main1( 3, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );		//14

	bb = 5.0;	Yleft = -2.5; Yright = 2.5; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x5_slayer3.msh"; c0 = 0.5; C0 = 263.0/240.0; eps = 1e-12;							//INPUT
	main1( 3, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );		//10

	bb = 1.0;	Yleft = -0.5; Yright = 0.5; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x1_slayer3.msh"; c0 = 0.5; C0 = 0.0; eps = 1e-12;							//INPUT
	main1( 3, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );		//5
*/

/*
	bb = 1.0;	Yleft = -0.5; Yright = 0.5; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x1_slayer3.msh"; c0 = 0.5; C0 = 0.0; eps = 1e-12;							//INPUT
	main1( 0, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 1, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 2, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 3, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 4, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );		//5

	bb = 5.0;	Yleft = -2.5; Yright = 2.5; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x5_slayer3.msh"; c0 = 0.5; C0 = 263.0/240.0; eps = 1e-12;							//INPUT
	main1( 0, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 1, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 2, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 3, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 4, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );		//10

	bb = 10.0;	Yleft = -5.0; Yright = 5.0; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x10_slayer3.msh"; c0 = 0.5; C0 = 4.325; eps = 1e-12;							//INPUT
	main1( 0, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 1, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 2, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 3, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );		//14


	bb = 10.0;	Yleft = -5.0; Yright = 5.0; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x10_slayer3_unregular.msh"; c0 = 0.5; C0 = 4.325; eps = 1e-12;							//INPUT
	main1( 1, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 0, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );		//16

	bb = 5.0;	Yleft = -2.5; Yright = 2.5; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x5_slayer3_unregular.msh"; c0 = 0.5; C0 = 263.0/240.0; eps = 1e-12;							//INPUT
	main1( 1, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 0, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );		//18

	bb = 1.0;	Yleft = -0.5; Yright = 0.5; Xdown = 0.0; Xup = 1.0; FileOfGrid = "net/1x1_slayer3_unregular.msh"; c0 = 0.5; C0 = 0.0; eps = 1e-12;							//INPUT
	main1( 1, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );
	main1( 0, bb, Yleft, Yright, Xdown, Xup, FileOfGrid, Yung1, Yung2, Yung3, c0, C0, eps );		//20
*/
}


int main()
{
	double Yung1		= 0.0;
	double Yung2		= 0.0;
	double Yung3		= 0.0;


/*
	Yung1		= 12.2;							//INPUT
	Yung2		= 1.2;							//INPUT
	Yung3		= 12.2;							//INPUT
	main2(Yung1, Yung2, Yung3);

	Yung1		= 8.8;							//INPUT
	Yung2		= 1.0;							//INPUT
	Yung3		= 10.2;							//INPUT
	main2(Yung1, Yung2, Yung3);

 	Yung1		= 3.3;							//INPUT
	Yung2		= 12.3;							//INPUT
	Yung3		= 6.1;							//INPUT
	main2(Yung1, Yung2, Yung3);

	Yung1		= 7.8;							//INPUT
	Yung2		= 3.6;							//INPUT
	Yung3		= 1.2;							//INPUT
	main2(Yung1, Yung2, Yung3);
*/

	Yung1		= 1.0;							//INPUT
	Yung2		= 1.0;							//INPUT
	Yung3		= 1.0;							//INPUT
	main2(Yung1, Yung2, Yung3);



}





