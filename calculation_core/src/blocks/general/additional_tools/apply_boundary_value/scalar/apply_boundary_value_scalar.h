#ifndef apply_boundary_value_scalar_def
#define apply_boundary_value_scalar_def 1

#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include "../../../domain/domain.h"
#include "../../../system_linear_algebraic_equations/system_linear_algebraic_equations.h"
#include "../../../boundary_value/boundary_value.h"

//! Дополнительный инструментарий
/*!
 * Additional tools \n
 * Сюда входят инструменты (функции) не вошедшие в другие, специализированные
 * пространства.
 */
namespace ATools
{
    //! Добавление к СЛАУ граничного скалярного условия
    /*!
     * Выполняется в таком виде: \n
     * ATools ::apply_boundary_value_scalar (boundary_value) .to_slae (slae, domain);
    */
    template<u8 dim>
    class apply_boundary_value_scalar
    {
        public:
            apply_boundary_value_scalar (const BoundaryValueScalar<dim> &bv) : 
                boundary_value (bv) {};

            void to_slae (SystemsLinearAlgebraicEquations &slae, 
                          const Domain<dim> &domain) const
            {
                if (boundary_value.boundary_type _is TBV::Dirichlet)
                {
                    std::map<u32, dbl> list_boundary_values;

                    dealii::VectorTools::interpolate_boundary_values (
                            domain.dof_handler,
                            boundary_value.boundary_id,
                            boundary_value.function,
                            list_boundary_values);

                    dealii::MatrixTools::apply_boundary_values (
                            list_boundary_values,
                            slae.matrix,
                            slae.solution,
                            slae.rhsv);
                }
                else if (boundary_value.boundary_type _is TBV::Neumann)
                {
                    dealii::Vector<dbl> tmp (slae.rhsv.size());

                    std::set<u8> b_id;
                    b_id.insert(boundary_value.boundary_id);

                    dealii::VectorTools::create_boundary_right_hand_side (
                            domain.dof_handler,
                            dealii::QGauss<dim-1>(2),
                            boundary_value.function,
                            tmp,
                            b_id);
                    // dealii::VectorTools::create_boundary_right_hand_side (
                    //         domain.dof_handler,
                    //         dealii::QGauss<dim-1>(2),
                    //         boundary_value.function,
                    //         tmp);

                    slae.rhsv += tmp;
                    // slae.rhsv = tmp;
                };
            };
        private:
            apply_boundary_value_scalar () {};
            BoundaryValueScalar<dim> boundary_value;
    };
};

#endif
