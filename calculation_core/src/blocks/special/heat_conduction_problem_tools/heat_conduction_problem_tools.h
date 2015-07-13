#ifndef heat_conduction_problem_def
#define heat_conduction_problem_def 1

#include <fstream>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/numerics/data_out.h>

#include "../../general/4_points_function/4_points_function.h"

//! Инструменты для решения задачи теплопроводности
namespace HCPTools
{
    //! Задать тензор теплопроводности 
    template<u8 dim>
    void set_thermal_conductivity(
            arr<arr<vec<dbl>, dim>, dim> &tensor, const arr<arr<vec<dbl>, dim>, dim> &coef)
    {
        for (st i = 0; i < dim; ++i)
        {
            for (st j = 0; j < dim; ++j)
            {
                tensor[i][j] .resize (coef[i][j].size());
                for (st k = 0; k < coef[i][j].size(); ++k)
                {
                    tensor[i][j][k] = coef[i][j][k];
                };
            };
        };
    };

    //! Распечатать в файл температуру
    template<u8 dim>
    void print_temperature (const dealii::Vector<dbl> &temperature, 
                            const dealii::DoFHandler<dim> &dof_handler,
                            const str file_name,
                            const dealii::DataOutBase::OutputFormat output_format = dealii::DataOutBase::gnuplot)
    {
        dealii::DataOut<dim> data_out;
        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (temperature, "temperature");
        data_out.build_patches ();

        auto name = file_name;
        // name += ".gpd";

        std::ofstream output (name);
        data_out.write (output, output_format);
    };

    //! Распечатать в файл 2d срез температуры для 3d случая
    void print_temperature_slice (const dealii::Vector<dbl> &temperature, 
                            const dealii::DoFHandler<3> &dof_handler,
                            const str file_name,
                            cst ort,
                            cdbl slice_coor)
    {
        FILE* f_out;
        f_out = fopen (file_name.c_str(), "w");
        for (auto cell = dof_handler.begin_active (); cell != dof_handler.end (); ++cell)
        {
            for (st i = 0; i < dealii::GeometryInfo<3>::vertices_per_cell; ++i)
            {
                if (std::abs(cell->vertex(i)(ort) - slice_coor) < 1e-10)
                {
                    fprintf(f_out, "%f %f %f %f\n",
                            cell->vertex(i)(0),
                            cell->vertex(i)(1),
                            cell->vertex(i)(2),
                            temperature(cell->vertex_dof_index(i, 0)));
                };
            };
        };
        fclose(f_out);
    };

    //! Распечатать в файл тепловые потоки
    template<u8 dim>
    void print_heat_conductions (const dealii::Vector<dbl> &temperature, 
                                 vec<arr<arr<dbl, 3>, 3>> &coef,
                                 const Domain<2> &domain,
                                 const str file_name)
    {
    };

    template<>
    void print_heat_conductions<2> (const dealii::Vector<dbl> &temperature, 
                                    vec<arr<arr<dbl, 3>, 3>> &coef,
                                    const Domain<2> &domain,
                                    const str file_name)
    {
        cu32 dofs_per_cell = domain.dof_handler.get_fe().dofs_per_cell;

        arr<dealii::Vector<dbl>, 2> grad;
        arr<dealii::Vector<dbl>, 2> hc;
        
        grad[0] .reinit (domain.dof_handler.n_dofs());
        grad[1] .reinit (domain.dof_handler.n_dofs());

        hc[0] .reinit (domain.dof_handler.n_dofs());
        hc[1] .reinit (domain.dof_handler.n_dofs());

        vec<st> N(domain.dof_handler.n_dofs());

        for (
                auto cell = domain.dof_handler.begin_active(); 
                cell     != domain.dof_handler.end(); 
                ++cell
            )
        {
            /* Точки 3 и 2 переставленны местами, потому что в диле у них
             * порядок зигзагом, а мне надо по кругу
            */
            arr<prmt::Point<2>, 4> points = {
                prmt::Point<2>(cell->vertex(0)(0), cell->vertex(0)(1)),
                prmt::Point<2>(cell->vertex(1)(0), cell->vertex(1)(1)),
                prmt::Point<2>(cell->vertex(3)(0), cell->vertex(3)(1)),
                prmt::Point<2>(cell->vertex(2)(0), cell->vertex(2)(1))};
            arr<dbl, 4> values = {
                temperature(cell->vertex_dof_index (0, 0)),
                temperature(cell->vertex_dof_index (1, 0)),
                temperature(cell->vertex_dof_index (3, 0)),
                temperature(cell->vertex_dof_index (2, 0))};

            Scalar4PointsFunc<2> function_on_cell(points, values);

            for (st i = 0; i < dofs_per_cell; ++i)
            {
                auto indx = cell->vertex_dof_index(i, 0);

                grad[0][indx] += function_on_cell.dx(cell->vertex(i));
                grad[1][indx] += function_on_cell.dy(cell->vertex(i));
                ++(N[indx]); 

                hc[0][indx] = 
                    coef[cell->material_id()][0][0] * grad[0][indx] +
                    coef[cell->material_id()][0][1] * grad[1][indx];
                hc[1][indx] = 
                    coef[cell->material_id()][1][0] * grad[0][indx] +
                    coef[cell->material_id()][1][1] * grad[1][indx];
            };
        };
        for (st i = 0; i < N.size(); ++i)
        {
            hc[0][i] /= N[i];
            hc[1][i] /= N[i];
        };

        {
            dealii::DataOut<2> data_out;
            data_out.attach_dof_handler (domain.dof_handler);
            data_out.add_data_vector (hc[0], "hc_dx");
            data_out.add_data_vector (hc[1], "hc_dy");
            data_out.build_patches ();

            auto name = file_name;

            std::ofstream output (name);
            data_out.write_gnuplot (output);
        };
    }

    //! Распечатать в файл тепловые градиенты
    template<u8 dim>
    void print_heat_gradient (const dealii::Vector<dbl> &temperature, 
                                 const Domain<2> &domain,
                                 const str file_name)
    {
    };

    template<>
    void print_heat_gradient<2> (const dealii::Vector<dbl> &temperature, 
                                    const Domain<2> &domain,
                                    const str file_name)
    {
        cu32 dofs_per_cell = domain.dof_handler.get_fe().dofs_per_cell;

        arr<dealii::Vector<dbl>, 2> grad;
        
        grad[0] .reinit (domain.dof_handler.n_dofs());
        grad[1] .reinit (domain.dof_handler.n_dofs());

        vec<st> N(domain.dof_handler.n_dofs());

        for (
                auto cell = domain.dof_handler.begin_active(); 
                cell     != domain.dof_handler.end(); 
                ++cell
            )
        {
            /* Точки 3 и 2 переставленны местами, потому что в диле у них
             * порядок зигзагом, а мне надо по кругу
            */
            arr<prmt::Point<2>, 4> points = {
                prmt::Point<2>(cell->vertex(0)(0), cell->vertex(0)(1)),
                prmt::Point<2>(cell->vertex(1)(0), cell->vertex(1)(1)),
                prmt::Point<2>(cell->vertex(3)(0), cell->vertex(3)(1)),
                prmt::Point<2>(cell->vertex(2)(0), cell->vertex(2)(1))};
            arr<dbl, 4> values = {
                temperature(cell->vertex_dof_index (0, 0)),
                temperature(cell->vertex_dof_index (1, 0)),
                temperature(cell->vertex_dof_index (3, 0)),
                temperature(cell->vertex_dof_index (2, 0))};

            Scalar4PointsFunc<2> function_on_cell(points, values);


            for (st i = 0; i < dofs_per_cell; ++i)
            {
                auto indx = cell->vertex_dof_index(i, 0);

                grad[0][indx] += function_on_cell.dx(cell->vertex(i));
                grad[1][indx] += function_on_cell.dy(cell->vertex(i));
                ++(N[indx]); 

            };
        };

        for (st i = 0; i < N.size(); ++i)
        {
            grad[0][i] /= N[i];
            grad[1][i] /= N[i];
        };

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//Вывод в файл
			int ID=-100;
			dealii::Vector<dbl> x_coordinate_of_vertex(domain.dof_handler.n_dofs());
			dealii::Vector<dbl> y_coordinate_of_vertex(domain.dof_handler.n_dofs());
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
					ID = cell->vertex_dof_index (1, 0);
					x_coordinate_of_vertex(ID) = cell->vertex(1)[0];
					y_coordinate_of_vertex(ID) = cell->vertex(1)[1];
						ID = cell->vertex_dof_index (2, 0);
						x_coordinate_of_vertex(ID) = cell->vertex(2)[0];
						y_coordinate_of_vertex(ID) = cell->vertex(2)[1];
							ID = cell->vertex_dof_index (3, 0);
							x_coordinate_of_vertex(ID) = cell->vertex(3)[0];
							y_coordinate_of_vertex(ID) = cell->vertex(3)[1];
			}
			FILE *TestFile1;
//			TestFile1 = fopen ("out/TestFile1.gpl", "w");
			TestFile1 = fopen ("out/GRAD.gpl", "w");
			for(int i=0; i<domain.dof_handler.n_dofs(); ++i)
			{
//																					№			"X"							"Y"						"F"			"dX"		 "dY"
				fprintf( TestFile1, "%d\t%5.35f\t%5.35f\t%5.35f\t%5.35f\t%5.35f\r\n",	i, x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), temperature(i), grad[0][i], grad[1][i] );
			}
			fclose(TestFile1);
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

/*
        {
            dealii::DataOut<2> data_out;
            data_out.attach_dof_handler (domain.dof_handler);
            data_out.add_data_vector (grad[0], "grad_dx");
            data_out.add_data_vector (grad[1], "drad_dy");
            data_out.build_patches ();

            auto name = file_name;

            std::ofstream output (name);
            data_out.write_gnuplot (output);
        };
*/


    }



    //! Распечатать в файл градиенты упругости
    template<u8 dim>
    void print_elasticity_gradient (const dealii::Vector<dbl> &elasticity, 
								const Domain<2> &domain,
								const str file_name)
    {
    };

    template<>
    void print_elasticity_gradient<2> (const dealii::Vector<dbl> &elasticity, 
								const Domain<2> &domain,
								const str file_name)
    {
        cu32 dofs_per_cell = domain.dof_handler.get_fe().dofs_per_cell;

        arr<dealii::Vector<dbl>, 2> grad;
        
        grad[0] .reinit (domain.dof_handler.n_dofs());
        grad[1] .reinit (domain.dof_handler.n_dofs());

        vec<st> N(domain.dof_handler.n_dofs());

        for (
                auto cell = domain.dof_handler.begin_active(); 
                cell     != domain.dof_handler.end(); 
                ++cell
            )
        {
            /* Точки 3 и 2 переставленны местами, потому что в диле у них
             * порядок зигзагом, а мне надо по кругу
            */
//------------------------------------------------------------------------------
//Вычисление произовдных для U_x
            arr<prmt::Point<2>, 4> points = {
                prmt::Point<2>(cell->vertex(0)(0), cell->vertex(0)(1)),			//координаты локальной вершины №0
                prmt::Point<2>(cell->vertex(1)(0), cell->vertex(1)(1)),			//координаты локальной вершины №1
                prmt::Point<2>(cell->vertex(3)(0), cell->vertex(3)(1)),			//координаты локальной вершины №2
                prmt::Point<2>(cell->vertex(2)(0), cell->vertex(2)(1))};		//координаты локальной вершины №3
            arr<dbl, 4> values = {
                elasticity(cell->vertex_dof_index(0, 0)),						//возвращет эл-нт вектора elasticity для вершины №0 для первой функции U_x
                elasticity(cell->vertex_dof_index(1, 0)),						//возвращет эл-нт вектора elasticity для вершины №1 для первой функции U_x
                elasticity(cell->vertex_dof_index(3, 0)),						//возвращет эл-нт вектора elasticity для вершины №2 для первой функции U_x
                elasticity(cell->vertex_dof_index(2, 0))};						//возвращет эл-нт вектора elasticity для вершины №3 для первой функции U_x

            Scalar4PointsFunc<2> function1_on_cell(points, values);
            for (st i = 0; i < dofs_per_cell; ++i)
            {
                auto indx = cell->vertex_dof_index(i, 0);						//возвращает глобальный DoF-индекс для i-ой вершины. 0 - первая функция U_x
                grad[0][indx] += function1_on_cell.dx(cell->vertex(i));			//передача координат x, y в функцию dx()
                grad[1][indx] += function1_on_cell.dy(cell->vertex(i));
                ++(N[indx]);
            };
//------------------------------------------------------------------------------
//Вычисление произовдных для U_y
//массив points - тот же самый, переписывать не станем
			values = {
                elasticity(cell->vertex_dof_index(0, 1)),						//возвращет эл-нт вектора elasticity для вершины №0 для первой функции U_y
                elasticity(cell->vertex_dof_index(1, 1)),						//возвращет эл-нт вектора elasticity для вершины №1 для первой функции U_y
                elasticity(cell->vertex_dof_index(3, 1)),						//возвращет эл-нт вектора elasticity для вершины №2 для первой функции U_y
                elasticity(cell->vertex_dof_index(2, 1))};						//возвращет эл-нт вектора elasticity для вершины №3 для первой функции U_y
            Scalar4PointsFunc<2> function2_on_cell(points, values);
            for (st i = 0; i < dofs_per_cell; ++i)
            {
                auto indx = cell->vertex_dof_index(i, 1);						//возвращает глобальный DoF-индекс для i-ой вершины. 0 - первая функция U_y
                grad[0][indx] += function2_on_cell.dx(cell->vertex(i));			//передача координат x, y в функцию dx()
                grad[1][indx] += function2_on_cell.dy(cell->vertex(i));
                ++(N[indx]);
            };

        };
//------------------------------------------------------------------------------
//поиск среднего арифметического значения производной в каждой вершине
        for (st i = 0; i < N.size(); ++i)
        {
            grad[0][i] /= N[i];
            grad[1][i] /= N[i];
//			std::cout << "grad\n";
//			std::cout << "grad[0][" << i << "] = " << grad[0][i] << "\n";
        };



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//Вывод в файл
			int ID=-100;
			dealii::Vector<dbl> x_coordinate_of_vertex(domain.dof_handler.n_dofs());
			dealii::Vector<dbl> y_coordinate_of_vertex(domain.dof_handler.n_dofs());
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
				ID = cell->vertex_dof_index (0, 1);
				x_coordinate_of_vertex(ID) = cell->vertex(0)[0];
				y_coordinate_of_vertex(ID) = cell->vertex(0)[1];
					ID = cell->vertex_dof_index (1, 0);
					x_coordinate_of_vertex(ID) = cell->vertex(1)[0];
					y_coordinate_of_vertex(ID) = cell->vertex(1)[1];
					ID = cell->vertex_dof_index (1, 1);
					x_coordinate_of_vertex(ID) = cell->vertex(1)[0];
					y_coordinate_of_vertex(ID) = cell->vertex(1)[1];
						ID = cell->vertex_dof_index (2, 0);
						x_coordinate_of_vertex(ID) = cell->vertex(2)[0];
						y_coordinate_of_vertex(ID) = cell->vertex(2)[1];
						ID = cell->vertex_dof_index (2, 1);
						x_coordinate_of_vertex(ID) = cell->vertex(2)[0];
						y_coordinate_of_vertex(ID) = cell->vertex(2)[1];
							ID = cell->vertex_dof_index (3, 0);
							x_coordinate_of_vertex(ID) = cell->vertex(3)[0];
							y_coordinate_of_vertex(ID) = cell->vertex(3)[1];
							ID = cell->vertex_dof_index (3, 1);
							x_coordinate_of_vertex(ID) = cell->vertex(3)[0];
							y_coordinate_of_vertex(ID) = cell->vertex(3)[1];
			}
			FILE *TestFile2;
//			TestFile2 = fopen ("out/TestFile2.gpl", "w");
			TestFile2 = fopen ("out/GRAD.gpl", "w");
			for(int i=0; i<domain.dof_handler.n_dofs(); ++i)
			{
//																					№			"X"							"Y"						"F"			"dX"		 "dY"
				fprintf( TestFile2, "%d\t%5.35f\t%5.35f\t%5.35f\t%5.35f\t%5.35f\r\n",	i, x_coordinate_of_vertex(i), y_coordinate_of_vertex(i), elasticity(i), grad[0][i], grad[1][i] );
			}
			fclose(TestFile2);
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
/*
        {
            dealii::DataOut<2> data_out;
            data_out.attach_dof_handler (domain.dof_handler);
            data_out.add_data_vector (grad[0], "grad_dx");
            data_out.add_data_vector (grad[1], "drad_dy");
            data_out.build_patches ();

            auto name = file_name;

            std::ofstream output (name);
            data_out.write_gnuplot (output);
        };
*/
    }







    //! Распечатать в файл тепловые градиенты
    template<u8 dim>
    double print_function_temperature (const dealii::Vector<dbl> &temperature, 
                                 const Domain<2> &domain,
								 dealii::Point<2, double> point)
    {
    };

    template<>
    double print_function_temperature <2> (const dealii::Vector<dbl> &temperature, 
                                    const Domain<2> &domain,
									dealii::Point<2, double> point)
    {
        cu32 dofs_per_cell = domain.dof_handler.get_fe().dofs_per_cell;

        arr<dealii::Vector<dbl>, 2> grad;

		int ii=0;
        for (
                auto cell = domain.dof_handler.begin_active(); 
                cell     != domain.dof_handler.end(); 
                ++cell
            )
        {
            /* Точки 3 и 2 переставленны местами, потому что в диле у них
             * порядок зигзагом, а мне надо по кругу
            */
            arr<prmt::Point<2>, 4> points = {
                prmt::Point<2>(cell->vertex(0)(0), cell->vertex(0)(1)),
                prmt::Point<2>(cell->vertex(1)(0), cell->vertex(1)(1)),
                prmt::Point<2>(cell->vertex(3)(0), cell->vertex(3)(1)),
                prmt::Point<2>(cell->vertex(2)(0), cell->vertex(2)(1))};
            arr<dbl, 4> values = {												//делим на 2 из-за того, что domain индексирует 2 функции поочерёдно
                temperature(cell->vertex_dof_index (0, 0)),
                temperature(cell->vertex_dof_index (1, 0)),
                temperature(cell->vertex_dof_index (3, 0)),
                temperature(cell->vertex_dof_index (2, 0))};
/*
			std::cout << "cell->vertex_dof_index (0, 0) / 2 = " << cell->vertex_dof_index (0, 0) / 2 << "\n";
			std::cout << "cell->vertex_dof_index (1, 0) / 2 = " << cell->vertex_dof_index (1, 0) / 2 << "\n";
			std::cout << "cell->vertex_dof_index (2, 0) / 2 = " << cell->vertex_dof_index (2, 0) / 2 << "\n";
			std::cout << "cell->vertex_dof_index (3, 0) / 2 = " << cell->vertex_dof_index (3, 0) / 2 << "\n";
			std::cout << "\n";
*/
			//проверка принадлежности переданной точки прямоугольнику:
			if(  ( (point(0)-points[0].x())*(point(0)-points[2].x()) <= 0.0 )&&( (point(1)-points[0].y())*(point(1)-points[2].y()) <= 0.0 )  )
			{
	            Scalar4PointsFunc<2> function_on_cell(points, values);
//				std::cout << "function_on_cell(point) = " << function_on_cell(point) << "\n";
//				std::cout << "ii = " << ii << "\n";
				return function_on_cell(point);
				break;
			}
			++ii;

        };
		std::cout << "ii = " << ii << "\n";
		
		return 0.0;
    }






    //! Распечатать в файл тепловые градиенты
    template<u8 dim>
    double print_function_elastic (const dealii::Vector<dbl> &Vec, 
                                 const Domain<2> &domain,
								 dealii::Point<2, double> point)
    {
    };

    template<>
    double print_function_elastic <2> (const dealii::Vector<dbl> &Vec, 
                                    const Domain<2> &domain,
									dealii::Point<2, double> point)
    {
        cu32 dofs_per_cell = domain.dof_handler.get_fe().dofs_per_cell;

        arr<dealii::Vector<dbl>, 2> grad;

		int ii=0;
        for (
                auto cell = domain.dof_handler.begin_active(); 
                cell     != domain.dof_handler.end(); 
                ++cell
            )
        {
            /* Точки 3 и 2 переставленны местами, потому что в диле у них
             * порядок зигзагом, а мне надо по кругу
            */
            arr<prmt::Point<2>, 4> points = {
                prmt::Point<2>(cell->vertex(0)(0), cell->vertex(0)(1)),
                prmt::Point<2>(cell->vertex(1)(0), cell->vertex(1)(1)),
                prmt::Point<2>(cell->vertex(3)(0), cell->vertex(3)(1)),
                prmt::Point<2>(cell->vertex(2)(0), cell->vertex(2)(1))};
            arr<dbl, 4> values = {												//делим на 2 из-за того, что domain индексирует 2 функции поочерёдно
                Vec(cell->vertex_dof_index (0, 0) / 2),
                Vec(cell->vertex_dof_index (1, 0) / 2),
                Vec(cell->vertex_dof_index (3, 0) / 2),
                Vec(cell->vertex_dof_index (2, 0) / 2)};
/*
			std::cout << "cell->vertex_dof_index (0, 0) / 2 = " << cell->vertex_dof_index (0, 0) / 2 << "\n";
			std::cout << "cell->vertex_dof_index (1, 0) / 2 = " << cell->vertex_dof_index (1, 0) / 2 << "\n";
			std::cout << "cell->vertex_dof_index (2, 0) / 2 = " << cell->vertex_dof_index (2, 0) / 2 << "\n";
			std::cout << "cell->vertex_dof_index (3, 0) / 2 = " << cell->vertex_dof_index (3, 0) / 2 << "\n";
			std::cout << "\n";
*/
			//проверка принадлежности переданной точки прямоугольнику:
			if(  ( (point(0)-points[0].x())*(point(0)-points[2].x()) <= 0.0 )&&( (point(1)-points[0].y())*(point(1)-points[2].y()) <= 0.0 )  )
			{
	            Scalar4PointsFunc<2> function_on_cell(points, values);
//				std::cout << "function_on_cell(point) = " << function_on_cell(point) << "\n";
//				std::cout << "ii = " << ii << "\n";
				return function_on_cell(point);
				break;
			}
			++ii;

        };
//		std::cout << "ii = " << ii << "\n";
		
		return 0.0;
    }





};

#endif
