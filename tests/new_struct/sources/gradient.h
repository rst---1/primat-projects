#ifndef GRADIENT
#define GRADIENT

#include <iostream>
#include <cmath>


#include <fstream>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/numerics/data_out.h>
//#include "../../general/4_points_function/4_points_function.h"

//создать файл GRAD.gpl - содержит производные
//			HCPTools::print_heat_gradient<2>(slae.solution, domain, "GRAD");
namespace GRADIENT
{

	using namespace std;


    void gradient (const dealii::Vector<dbl> &temperature, 
                                    const Domain<2> &domain)
    {

        cu32 dofs_per_cell = domain.dof_handler.get_fe().dofs_per_cell;

        
		dealii::Vector<dbl> NOV(domain.dof_handler.n_dofs());												//Вектор содержащий номера вершин. NOV - numbers of vertexes
		NOV = 0.0;
		dealii::Vector<dbl> x_coordinate_of_vertex(domain.dof_handler.n_dofs());
		dealii::Vector<dbl> y_coordinate_of_vertex(domain.dof_handler.n_dofs());
		dealii::Vector<dbl> gradX(domain.dof_handler.n_dofs());
		gradX = 0.0;
		dealii::Vector<dbl> gradY(domain.dof_handler.n_dofs());
		gradY = 0.0;
		




			int ID=-100;
            int indx = -100;
			double x_0, x_1, x_2, x_3, y_0, y_1, y_2, y_3;
			double f_0, f_1, f_2, f_3;
			double dx, dy, f_x, f_y;
			for (auto cell = domain.dof_handler.begin_active(); 
            cell != domain.dof_handler.end(); ++cell)
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
				x_coordinate_of_vertex(ID) = cell->vertex(0)[0];
				y_coordinate_of_vertex(ID) = cell->vertex(0)[1];
//				NOV[ID] = ID/2;
				NOV[ID] = ID;
				ID = cell->vertex_dof_index (0, 1);
				x_coordinate_of_vertex(ID) = cell->vertex(0)[0];
				y_coordinate_of_vertex(ID) = cell->vertex(0)[1];
//				NOV[ID] = (ID-1)/2;
				NOV[ID] = ID;
				x_0 = x_coordinate_of_vertex(ID); y_0 = y_coordinate_of_vertex(ID);
					ID = cell->vertex_dof_index (1, 0);
					x_coordinate_of_vertex(ID) = cell->vertex(1)[0];
					y_coordinate_of_vertex(ID) = cell->vertex(1)[1];
//					NOV[ID] = ID/2;
					NOV[ID] = ID;
					ID = cell->vertex_dof_index (1, 1);
					x_coordinate_of_vertex(ID) = cell->vertex(1)[0];
					y_coordinate_of_vertex(ID) = cell->vertex(1)[1];
//					NOV[ID] = (ID-1)/2;
					NOV[ID] = ID;
					x_1 = x_coordinate_of_vertex(ID); y_1 = y_coordinate_of_vertex(ID);
						ID = cell->vertex_dof_index (2, 0);
						x_coordinate_of_vertex(ID) = cell->vertex(2)[0];
						y_coordinate_of_vertex(ID) = cell->vertex(2)[1];
//						NOV[ID] = ID/2;
						NOV[ID] = ID;
						ID = cell->vertex_dof_index (2, 1);
						x_coordinate_of_vertex(ID) = cell->vertex(2)[0];
						y_coordinate_of_vertex(ID) = cell->vertex(2)[1];
//						NOV[ID] = (ID-1)/2;
						NOV[ID] = ID;
						x_2 = x_coordinate_of_vertex(ID); y_2 = y_coordinate_of_vertex(ID);
							ID = cell->vertex_dof_index (3, 0);
							x_coordinate_of_vertex(ID) = cell->vertex(3)[0];
							y_coordinate_of_vertex(ID) = cell->vertex(3)[1];
//							NOV[ID] = ID/2;
							NOV[ID] = ID;
							ID = cell->vertex_dof_index (3, 1);
							x_coordinate_of_vertex(ID) = cell->vertex(3)[0];
							y_coordinate_of_vertex(ID) = cell->vertex(3)[1];
//							NOV[ID] = (ID-1)/2;
							NOV[ID] = ID;
							x_3 = x_coordinate_of_vertex(ID); y_3 = y_coordinate_of_vertex(ID);
					
/*				std::cout << "0 - 1 :" << sqrt((x_0-x_1)*(x_0-x_1) + (y_0-y_1)*(y_0-y_1)) << "\n";
				std::cout << "2 - 3 :" << sqrt((x_2-x_3)*(x_2-x_3) + (y_2-y_3)*(y_2-y_3)) << "\n";
				std::cout << "2 - 0 :" << sqrt((x_2-x_0)*(x_2-x_0) + (y_2-y_0)*(y_2-y_0)) << "\n";
				std::cout << "3 - 1 :" << sqrt((x_3-x_1)*(x_3-x_1) + (y_3-y_1)*(y_3-y_1)) << "\n";
				std::cout << "0 - 3 :" << sqrt((x_0-x_3)*(x_0-x_3) + (y_0-y_3)*(y_0-y_3)) << "\n";
				std::cout << "2 - 1 :" << sqrt((x_2-x_1)*(x_2-x_1) + (y_2-y_1)*(y_2-y_1)) << "\n";
				std::cout << "\n";
*/
				f_0 = temperature(cell->vertex_dof_index (0, 0));				//вычисление значения функции U_x
				f_1 = temperature(cell->vertex_dof_index (1, 0));				//вычисление значения функции U_x
				f_2 = temperature(cell->vertex_dof_index (2, 0));				//вычисление значения функции U_x
				f_3 = temperature(cell->vertex_dof_index (3, 0));				//вычисление значения функции U_x
				if(x_0 == x_1)
				{
					dy = y_0 - y_1;												//вычисление длины грани dy
					dx = x_2 - x_0;												//вычисление длины грани dx

					f_x = (f_2 - f_0) / dx;										//вычисление производной
					ID = cell->vertex_dof_index (2, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;
					ID = cell->vertex_dof_index (0, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;

					f_x = (f_3 - f_1) / dx;										//вычисление производной
					ID = cell->vertex_dof_index (3, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;
					ID = cell->vertex_dof_index (1, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;

					f_y = (f_2 - f_3) / dy;										//вычисление производной
					ID = cell->vertex_dof_index (2, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;
					ID = cell->vertex_dof_index (3, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;											//вычисление производной

					f_y = (f_0 - f_1) / dy;										//вычисление производной
					ID = cell->vertex_dof_index (0, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;
					ID = cell->vertex_dof_index (1, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;											//вычисление производной
				}
				else
				{
					dx = x_0 - x_1;												//вычисление длины грани
					dy = y_2 - y_0;												//вычисление длины грани
					f_x = (f_0 - f_1) / dx;										//вычисление производной
					f_y = (f_2 - f_0) / dy;										//вычисление производной

					f_x = (f_3 - f_2) / dx;										//вычисление производной
					ID = cell->vertex_dof_index (3, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;
					ID = cell->vertex_dof_index (2, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;

					f_x = (f_1 - f_0) / dx;										//вычисление производной
					ID = cell->vertex_dof_index (1, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;
					ID = cell->vertex_dof_index (0, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;

					f_y = (f_1 - f_3) / dy;										//вычисление производной
					ID = cell->vertex_dof_index (1, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;
					ID = cell->vertex_dof_index (3, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;											//вычисление производной

					f_y = (f_0 - f_2) / dy;										//вычисление производной
					ID = cell->vertex_dof_index (0, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;
					ID = cell->vertex_dof_index (2, 0);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;											//вычисление производной
				}

				f_0 = temperature(cell->vertex_dof_index (0, 1));				//вычисление значения функции U_y
				f_1 = temperature(cell->vertex_dof_index (1, 1));				//вычисление значения функции U_y
				f_2 = temperature(cell->vertex_dof_index (2, 1));				//вычисление значения функции U_y
				f_3 = temperature(cell->vertex_dof_index (3, 1));				//вычисление значения функции U_y
				if(x_0 == x_1)
				{
					dy = y_0 - y_1;												//вычисление длины грани dy
					dx = x_2 - x_0;												//вычисление длины грани dx

					f_x = (f_2 - f_0) / dx;										//вычисление производной
					ID = cell->vertex_dof_index (2, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;
					ID = cell->vertex_dof_index (0, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;

					f_x = (f_3 - f_1) / dx;										//вычисление производной
					ID = cell->vertex_dof_index (3, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;
					ID = cell->vertex_dof_index (1, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;

					f_y = (f_2 - f_3) / dy;										//вычисление производной
					ID = cell->vertex_dof_index (2, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;
					ID = cell->vertex_dof_index (3, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;											//вычисление производной

					f_y = (f_0 - f_1) / dy;										//вычисление производной
					ID = cell->vertex_dof_index (0, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;
					ID = cell->vertex_dof_index (1, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;											//вычисление производной
				}
				else
				{
					dx = x_0 - x_1;												//вычисление длины грани
					dy = y_2 - y_0;												//вычисление длины грани
					f_x = (f_0 - f_1) / dx;										//вычисление производной
					f_y = (f_2 - f_0) / dy;										//вычисление производной

					f_x = (f_3 - f_2) / dx;										//вычисление производной
					ID = cell->vertex_dof_index (3, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;
					ID = cell->vertex_dof_index (2, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;

					f_x = (f_1 - f_0) / dx;										//вычисление производной
					ID = cell->vertex_dof_index (1, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;
					ID = cell->vertex_dof_index (0, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradX[ID] = f_x;

					f_y = (f_1 - f_3) / dy;										//вычисление производной
					ID = cell->vertex_dof_index (1, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;
					ID = cell->vertex_dof_index (3, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;											//вычисление производной

					f_y = (f_0 - f_2) / dy;										//вычисление производной
					ID = cell->vertex_dof_index (0, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;
					ID = cell->vertex_dof_index (2, 1);							//вычисление номера DoF, которому присвоим произвоодную
					gradY[ID] = f_y;											//вычисление производной
				}


			}

//------------------------------------------------------------------------------
//Вывод в файл
//			auto name = file_name;
//			name += ".gpd";
			FILE *outFile;
			outFile = fopen ("out/GRAD.gpl", "w");
			for(int i=0; i<domain.dof_handler.n_dofs(); ++i)
				{
//slae.solution(i) - это массив, с чередованием функций {U_x, U_y, U_x, U_y..}
//Поэтому далее, в некоторых случаях: 		U_x = slae.solution(i),
//									 		U_y = slae.solution(i+1)
//																							"x"						"y"				 "NOV"
//					fprintf(  outFile,"%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\r\n",	x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), NOV(i) );
					fprintf(  outFile,"%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), NOV(i), gradX(i), gradY(i) );
//					!(i % 2) ? fprintf(  outFile,"%9.4f\t%9.4f\t%9.4f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), tau_xx(i) )
//							 : fprintf(  outFile,"%9.4f\t%9.4f\t%9.4f\r\n", x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), tau_yy(i) );

				}
		fclose(outFile);
//------------------------------------------------------------------------------






/*
        FILE *F1;
        F1 = fopen("gradx.gpd", "w");
        for (st i = 0; i < ps.size(); ++i)
        {
            fprintf(F1, "%lf %lf %lf\n", ps[i].x(), ps[i].y(), gradx[i]);
        };
        fclose(F1);
*/
    }

}


#endif
