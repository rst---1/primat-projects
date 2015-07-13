#ifndef nikola_elastic_problem_def
#define nikola_elastic_problem_def 1

#include <fstream>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/numerics/data_out.h>

//#include "../../general/4_points_function/4_points_function.h"

//! Инструменты для решения задачи теплопроводности
namespace NEPTools
{
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

};

#endif
