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


#include "../../../calculation_core/src/blocks/special/nikola_problem/source/scalar/source_scalar.h"

#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/grid/grid_reordering.h>

#include <deal.II/grid/grid_tools.h>						//для gmsh
#include <deal.II/grid/grid_in.h>						//для gmsh

extern void make_grid(
        dealii::Triangulation< 2 >&,
        vec<prmt::Point<2>>,
        vec<st>);

extern void set_grid(
        dealii::Triangulation< 2 >&,
        vec<prmt::Point<2>>,
        vec<prmt::Point<2>>);

void debputs()
{
    static int n = 0;
    printf("DEBUG %d\n", n);
    n++;
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
class Beeline
{
    public:
        Beeline (arr<prmt::Point<2>, 4> nodes, arr<dbl, 4> f_values)
        {
            cdbl x1 = nodes[0].x();
            cdbl x2 = nodes[1].x();
            cdbl x3 = nodes[2].x();
            cdbl x4 = nodes[3].x();

            cdbl y1 = nodes[0].y();
            cdbl y2 = nodes[1].y();
            cdbl y3 = nodes[2].y();
            cdbl y4 = nodes[3].y();

            cdbl f1 = f_values[0]; 
            cdbl f2 = f_values[1];
            cdbl f3 = f_values[2];
            cdbl f4 = f_values[3];

            a=-(x2*y2*x3*y4*f1+y1*x3*x4*y4*f2-y1*x3*x2*y2*f4-x1*y1*x3*y4*f2-x1*y3*x4*y4*f2-x1*f3*x2*y2*y4+x1*x3*y3*f2*y4+x1*y3*x2*y2*f4+y1*f3*x4*x2*y2-y3*x4*x2*y2*f1-y3*y1*x3*x4*f2+x1*y3*x4*y1*f2-x2*y3*x3*y4*f1+x1*y1*f3*y4*x2+y1*y3*x3*x2*f4+x2*y3*x4*y4*f1-y1*f3*y4*x4*x2+x3*y3*x4*y2*f1-x3*f1*y2*x4*y4+y2*x3*x1*y1*f4+y2*x1*f3*y4*x4-y2*x4*x1*y1*f3-y2*x1*y3*x3*f4-x1*y1*y3*x2*f4)/(-y1*y3*x3*x2-x4*x2*y2*y1-x3*y3*x4*y2+x2*y3*x3*y4-x2*y2*x3*y4+x2*y1*x4*y4-x2*y3*x4*y4+y3*x4*x2*y2+y1*x3*x2*y2+y2*x3*x4*y4-y1*x3*x4*y4+x1*y1*x3*y4+y2*x1*y3*x3-x1*y4*x2*y1+x1*x2*y2*y4-y2*x3*x1*y1+y3*y1*x3*x4+x1*y3*x4*y4-x1*y3*x2*y2+x1*y1*y3*x2-x1*x3*y3*y4-x1*y3*x4*y1-x1*y2*x4*y4+y2*x4*x1*y1);
            b=-(x1*y1*y2*f3-x1*y1*y2*f4-x1*y1*f3*y4+x1*y1*y4*f2-x1*y1*f2*y3+x1*y1*y3*f4+y3*x3*y2*f4-y2*x3*y3*f1+y3*x4*y4*f2-y3*x4*y4*f1+x3*y3*f1*y4-x3*y3*f2*y4+f3*x2*y2*y4-x2*y2*f1*y4+x2*y2*f1*y3-y3*x2*y2*f4-y1*x4*y4*f2-y1*y3*x3*f4+y1*x2*y2*f4+y1*x3*y3*f2-y1*x2*y2*f3-f3*y2*x4*y4+y1*f3*y4*x4+f1*y2*x4*y4)/(-y1*y3*x3*x2-x4*x2*y2*y1-x3*y3*x4*y2+x2*y3*x3*y4-x2*y2*x3*y4+x2*y1*x4*y4-x2*y3*x4*y4+y3*x4*x2*y2+y1*x3*x2*y2+y2*x3*x4*y4-y1*x3*x4*y4+x1*y1*x3*y4+y2*x1*y3*x3-x1*y4*x2*y1+x1*x2*y2*y4-y2*x3*x1*y1+y3*y1*x3*x4+x1*y3*x4*y4-x1*y3*x2*y2+x1*y1*y3*x2-x1*x3*y3*y4-x1*y3*x4*y1-x1*y2*x4*y4+y2*x4*x1*y1);
            c=(x1*x2*y2*f4-x1*f3*x2*y2+x3*x1*y1*f4-x1*y1*x2*f4+x1*y1*x2*f3-x1*x4*y4*f2+x4*x1*y1*f2+x1*f3*y4*x4-x4*x1*y1*f3+x1*x3*y3*f2-x3*x1*y1*f2-x1*y3*x3*f4-x3*x2*y2*f4+x3*x2*y2*f1-x4*x2*y2*f1+x4*x2*y2*f3-f3*y4*x4*x2+x2*x4*y4*f1+y3*x3*x2*f4-x4*x3*y3*f2-x2*x3*y3*f1+x3*x4*y4*f2-x3*x4*y4*f1+x4*x3*y3*f1)/(-y1*y3*x3*x2-x4*x2*y2*y1-x3*y3*x4*y2+x2*y3*x3*y4-x2*y2*x3*y4+x2*y1*x4*y4-x2*y3*x4*y4+y3*x4*x2*y2+y1*x3*x2*y2+y2*x3*x4*y4-y1*x3*x4*y4+x1*y1*x3*y4+y2*x1*y3*x3-x1*y4*x2*y1+x1*x2*y2*y4-y2*x3*x1*y1+y3*y1*x3*x4+x1*y3*x4*y4-x1*y3*x2*y2+x1*y1*y3*x2-x1*x3*y3*y4-x1*y3*x4*y1-x1*y2*x4*y4+y2*x4*x1*y1);
            d=(-x3*y1*f4+x3*y1*f2+y1*f3*x4-x4*y1*f2-x1*y2*f4+x3*y2*f4-x3*f1*y2+y2*x1*f3+f1*y2*x4-f3*y2*x4+x4*y3*f2-x1*y3*f2-x4*y3*f1+x1*y4*f2-x3*y4*f2+x3*y4*f1-x2*f1*y4+f3*x2*y4+x2*f1*y3-y3*x2*f4+x1*y3*f4-x1*f3*y4+y1*x2*f4-y1*x2*f3)/(-y1*y3*x3*x2-x4*x2*y2*y1-x3*y3*x4*y2+x2*y3*x3*y4-x2*y2*x3*y4+x2*y1*x4*y4-x2*y3*x4*y4+y3*x4*x2*y2+y1*x3*x2*y2+y2*x3*x4*y4-y1*x3*x4*y4+x1*y1*x3*y4+y2*x1*y3*x3-x1*y4*x2*y1+x1*x2*y2*y4-y2*x3*x1*y1+y3*y1*x3*x4+x1*y3*x4*y4-x1*y3*x2*y2+x1*y1*y3*x2-x1*x3*y3*y4-x1*y3*x4*y1-x1*y2*x4*y4+y2*x4*x1*y1);
        };
        dbl operator () (cdbl x, cdbl y) 
        {
            return a + b * x + c * y + d * x * y;
        };
        dbl dx (cdbl x, cdbl y) 
        {
            return b + d * y;
        };
        dbl dy (cdbl x, cdbl y) 
        {
            return b + d * y;
        };
        dbl operator () (const prmt::Point<2> &p) 
        {
            return operator()(p.x(), p.y());
        };
        dbl dx (const prmt::Point<2> &p) 
        {
            return dx(p.x(), p.y());
        };
        dbl dy (const prmt::Point<2> &p)
        {
            return dy(p.x(), p.y());
        };

        dbl a, b, c, d;
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
template<uint8_t dim>
dealii::Point<dim, double> get_grad (
        const typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
        const dealii::Vector<double> solution,
        uint8_t index_vertex)
{
    dealii::Point<dim, double> grad;

    double x1 = cell->vertex(0)(0);
    double x2 = cell->vertex(1)(0);
    double x3 = cell->vertex(2)(0);
    double x4 = cell->vertex(3)(0);

    double y1 = cell->vertex(0)(1);
    double y2 = cell->vertex(1)(1);
    double y3 = cell->vertex(2)(1);
    double y4 = cell->vertex(3)(1);

    double f1 = solution(cell->vertex_dof_index (0, 0));
    double f2 = solution(cell->vertex_dof_index (1, 0));
    double f3 = solution(cell->vertex_dof_index (2, 0));
    double f4 = solution(cell->vertex_dof_index (3, 0));

    double b=-(x1*y1*y2*f3-x1*y1*y2*f4-x1*y1*f3*y4+x1*y1*y4*f2-x1*y1*f2*y3+x1*y1*y3*f4+y3*x3*y2*f4-y2*x3*y3*f1+y3*x4*y4*f2-y3*x4*y4*f1+x3*y3*f1*y4-x3*y3*f2*y4+f3*x2*y2*y4-x2*y2*f1*y4+x2*y2*f1*y3-y3*x2*y2*f4-y1*x4*y4*f2-y1*y3*x3*f4+y1*x2*y2*f4+y1*x3*y3*f2-y1*x2*y2*f3-f3*y2*x4*y4+y1*f3*y4*x4+f1*y2*x4*y4)/(-y1*y3*x3*x2-x4*x2*y2*y1-x3*y3*x4*y2+x2*y3*x3*y4-x2*y2*x3*y4+x2*y1*x4*y4-x2*y3*x4*y4+y3*x4*x2*y2+y1*x3*x2*y2+y2*x3*x4*y4-y1*x3*x4*y4+x1*y1*x3*y4+y2*x1*y3*x3-x1*y4*x2*y1+x1*x2*y2*y4-y2*x3*x1*y1+y3*y1*x3*x4+x1*y3*x4*y4-x1*y3*x2*y2+x1*y1*y3*x2-x1*x3*y3*y4-x1*y3*x4*y1-x1*y2*x4*y4+y2*x4*x1*y1);
    double c=(x1*x2*y2*f4-x1*f3*x2*y2+x3*x1*y1*f4-x1*y1*x2*f4+x1*y1*x2*f3-x1*x4*y4*f2+x4*x1*y1*f2+x1*f3*y4*x4-x4*x1*y1*f3+x1*x3*y3*f2-x3*x1*y1*f2-x1*y3*x3*f4-x3*x2*y2*f4+x3*x2*y2*f1-x4*x2*y2*f1+x4*x2*y2*f3-f3*y4*x4*x2+x2*x4*y4*f1+y3*x3*x2*f4-x4*x3*y3*f2-x2*x3*y3*f1+x3*x4*y4*f2-x3*x4*y4*f1+x4*x3*y3*f1)/(-y1*y3*x3*x2-x4*x2*y2*y1-x3*y3*x4*y2+x2*y3*x3*y4-x2*y2*x3*y4+x2*y1*x4*y4-x2*y3*x4*y4+y3*x4*x2*y2+y1*x3*x2*y2+y2*x3*x4*y4-y1*x3*x4*y4+x1*y1*x3*y4+y2*x1*y3*x3-x1*y4*x2*y1+x1*x2*y2*y4-y2*x3*x1*y1+y3*y1*x3*x4+x1*y3*x4*y4-x1*y3*x2*y2+x1*y1*y3*x2-x1*x3*y3*y4-x1*y3*x4*y1-x1*y2*x4*y4+y2*x4*x1*y1);
    double d=(-x3*y1*f4+x3*y1*f2+y1*f3*x4-x4*y1*f2-x1*y2*f4+x3*y2*f4-x3*f1*y2+y2*x1*f3+f1*y2*x4-f3*y2*x4+x4*y3*f2-x1*y3*f2-x4*y3*f1+x1*y4*f2-x3*y4*f2+x3*y4*f1-x2*f1*y4+f3*x2*y4+x2*f1*y3-y3*x2*f4+x1*y3*f4-x1*f3*y4+y1*x2*f4-y1*x2*f3)/(-y1*y3*x3*x2-x4*x2*y2*y1-x3*y3*x4*y2+x2*y3*x3*y4-x2*y2*x3*y4+x2*y1*x4*y4-x2*y3*x4*y4+y3*x4*x2*y2+y1*x3*x2*y2+y2*x3*x4*y4-y1*x3*x4*y4+x1*y1*x3*y4+y2*x1*y3*x3-x1*y4*x2*y1+x1*x2*y2*y4-y2*x3*x1*y1+y3*y1*x3*x4+x1*y3*x4*y4-x1*y3*x2*y2+x1*y1*y3*x2-x1*x3*y3*y4-x1*y3*x4*y1-x1*y2*x4*y4+y2*x4*x1*y1);

    grad(0) = b + d * cell->vertex(index_vertex)(1);
    grad(1) = c + d * cell->vertex(index_vertex)(0);

    return grad;
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
template <size_t num_points>
void set_tria(dealii::Triangulation< 2 > &triangulation, 
        const double points[num_points], 
        const size_t material_id[num_points - 1][num_points - 1])
{

    const size_t num_cells = num_points - 1;

    std::vector< dealii::Point< 2 > > v (num_points * num_points);

    FOR_I (0, num_points)
        FOR_J (0, num_points)
        {
            v[i * num_points + j] = dealii::Point< 2 >(points[j], points[i]);
        };

    std::vector< dealii::CellData< 2 > > c (
            num_cells * num_cells, dealii::CellData< 2 >());

    FOR_I (0, num_cells)
        FOR_J (0, num_cells)
        {
            c[i * num_cells + j].vertices[0] = i * num_points + j + 0;
            c[i * num_cells + j].vertices[1] = i * num_points + j + 1;
            c[i * num_cells + j].vertices[2] = i * num_points + j + num_points;
            c[i * num_cells + j].vertices[3] = i * num_points + j + num_points + 1;

            c[i * num_cells + j].material_id = material_id[i][j];
        };

    triangulation .create_triangulation (v, c, dealii::SubCellData());
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
template <u8 dim>
void solved_heat_problem_on_cell (
        const dealii::Triangulation<dim> &grid,
        const vec<ATools::SecondOrderTensor> &coef,
        OnCell::SystemsLinearAlgebraicEquations<dim> &slae)
{
    enum {x, y, z};

    Domain<dim> domain;
    domain.grid .copy_triangulation(grid);
    dealii::FE_Q<dim> fe(1);
    domain.dof_init (fe);

    OnCell::BlackOnWhiteSubstituter bows;

    LaplacianScalar<dim> element_matrix (domain.dof_handler.get_fe());

    element_matrix.C .resize(2);
    for (st i = 0; i < coef.size(); ++i)
    {
        element_matrix.C[i][x][x] = coef[i][x][x];
        element_matrix.C[i][x][y] = coef[i][x][y];
        element_matrix.C[i][y][x] = coef[i][y][x];
        element_matrix.C[i][y][y] = coef[i][y][y];
    };

    const bool scalar_type = 0;
    OnCell::prepare_system_equations<scalar_type> (slae, bows, domain);

    OnCell::Assembler::assemble_matrix<dim> (slae.matrix, element_matrix, domain.dof_handler, bows);

    FOR(i, 0, dim)
    {
        vec<arr<dbl, 2>> coef_for_rhs(2);
        FOR(j, 0, element_matrix.C.size())
        {
            FOR(k, 0, 2)
            {
                coef_for_rhs[j][k] = element_matrix.C[j][i][k];
            };
        };
        OnCell::SourceScalar<dim> element_rhsv (coef_for_rhs, domain.dof_handler.get_fe());
        OnCell::Assembler::assemble_rhsv<dim> (slae.rhsv[i], element_rhsv, domain.dof_handler, bows);
        {
            dealii::DataOut<dim> data_out;
            data_out.attach_dof_handler (domain.dof_handler);
            data_out.add_data_vector (slae.rhsv[0], "xb");
            data_out.add_data_vector (slae.rhsv[1], "yb");
            data_out.build_patches ();

            auto name = "b.gpd";

            std::ofstream output (name);
            data_out.write_gnuplot (output);
        };

        dealii::SolverControl solver_control (10000, 1e-12);
        dealii::SolverCG<> solver (solver_control);
        solver.solve (
                slae.matrix,
                slae.solution[i],
                slae.rhsv[i]
                ,dealii::PreconditionIdentity()
                );
        FOR(j, 0, slae.solution[i].size())
            slae.solution[i][j] = slae.solution[i][bows.subst (j)];
    };
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
template<size_t size>
std::array<size_t, 2> to2D(const size_t i)//, const size_t j)
{
    std::array<size_t, 2> res;

    switch (size)
    {
        case 6*6:
            {
                if (i < 3)
                {
                    res[0] = i;
                    res[1] = i;
                }
                else
                {
                    switch (i)
                    {
                        case 3: res[0]=1; res[1]=2; break; 
                        case 4: res[0]=2; res[1]=0; break; 
                        case 5: res[0]=0; res[1]=1; break;
                    };
                };
            };
            break;
        case 9*9:
            {
                    res[0] = i / 3;
                    res[1] = i % 3;
            };
            break;
    };

    return res;
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
template<size_t size>
void print_tensor(const ATools::FourthOrderTensor &tensor)
{
    const size_t width = static_cast<size_t>(sqrt(size));

    for (size_t i = 0; i < width; ++i)
    {
        auto ind = to2D<size>(i);
        uint8_t im = ind[0];
        uint8_t in = ind[1];

        for (size_t j = 0; j < width; ++j)
        {
            auto jnd = to2D<size>(j);
            uint8_t jm = jnd[0];
            uint8_t jn = jnd[1];

            if (fabs(tensor[im][in][jm][jn]) > 0.0000001)
                printf("\x1B[31m%f\x1B[0m   ", 
                        tensor[im][in][jm][jn]);
            else
                printf("%f   ", 
                        tensor[im][in][jm][jn]);
        };
        for (size_t i = 0; i < 2; ++i)
            printf("\n");
    };

    printf("\n");
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
void set_hexagon_grid_pure(dealii::Triangulation< 2 > &triangulation, 
        const double len_edge,
        const double radius)
{
    double Ro = len_edge;
    double ro = radius;
    double Ri = Ro * (sqrt(3.0) / 2.0);
    double ri = ro * (sqrt(3.0) / 2.0);

    printf("Radius %f %f\n", ro, ri);

    double a[7] = {0.0, Ri-ri, Ri, Ri+ri, 2.0*Ri, ri, 2.0*Ri - ri};
    double b[15] = {
        0.0, ro/2.0, ro, Ro/2, Ro, 1.5*Ro - ro, 1.5*Ro - ro / 2.0,
        1.5*Ro, 1.5*Ro + ro / 2.0, 1.5*Ro + ro, 2.0*Ro, 2.5*Ro, 3.0*Ro-ro,
        3.0*Ro-ro/2.0, 3.0*Ro};


    std::vector<dealii::Point< 2 > > v (30); //30

    

    v[0][0]  = a[0]; v[0][1]  = b[0];
    v[1][0]  = a[1]; v[1][1]  = b[0];
    v[2][0]  = a[2]; v[2][1]  = b[0];
    v[3][0]  = a[3]; v[3][1]  = b[0];
    v[4][0]  = a[4]; v[4][1]  = b[0];

    v[5][0]  = a[1]; v[5][1]  = b[1];
    v[6][0]  = a[3]; v[6][1]  = b[1];

    v[7][0]  = a[2]; v[7][1]  = b[2];

    v[8][0]  = a[0]; v[8][1]  = b[3];
    v[9][0]  = a[4]; v[9][1]  = b[3];

    v[10][0] = a[2]; v[10][1] = b[4];

    v[11][0] = a[0]; v[11][1] = b[5];
    v[12][0] = a[4]; v[12][1] = b[5];

    v[13][0] = a[5]; v[13][1] = b[6];
    v[14][0] = a[6]; v[14][1] = b[6];

    v[15][0] = a[5]; v[15][1] = b[8];
    v[16][0] = a[6]; v[16][1] = b[8];

    v[17][0] = a[0]; v[17][1] = b[9];
    v[18][0] = a[4]; v[18][1] = b[9];

    v[19][0] = a[2]; v[19][1] = b[10];

    v[20][0] = a[0]; v[20][1] = b[11];
    v[21][0] = a[4]; v[21][1] = b[11];

//    v[23][0] = a[4]; v[23][1] = b[8];

    v[22][0] = a[2]; v[22][1] = b[12];

    v[23][0] = a[1]; v[23][1] = b[13];
    v[24][0] = a[3]; v[24][1] = b[13];

    v[25][0] = a[0]; v[25][1] = b[14];
    v[26][0] = a[1]; v[26][1] = b[14];
    v[27][0] = a[2]; v[27][1] = b[14];
    v[28][0] = a[3]; v[28][1] = b[14];
    v[29][0] = a[4]; v[29][1] = b[14];
//
////    v[31][0] = a[4]; v[31][1] = b[8];  // 13

    std::vector< dealii::CellData< 2 > > c (20, dealii::CellData<2>()); //20


    c[0].vertices[0] = 1;
    c[0].vertices[1] = 5;
    c[0].vertices[2] = 8;
    c[0].vertices[3] = 0;
    c[0].material_id = 0;

    c[1].vertices[0] = 1;
    c[1].vertices[1] = 2;
    c[1].vertices[2] = 7;
    c[1].vertices[3] = 5;
    c[1].material_id = 1;

    c[2].vertices[0] = 2;
    c[2].vertices[1] = 3;
    c[2].vertices[2] = 6;
    c[2].vertices[3] = 7;
    c[2].material_id = 1;

    c[3].vertices[0] = 3;
    c[3].vertices[1] = 4;
    c[3].vertices[2] = 9;
    c[3].vertices[3] = 6;
    c[3].material_id = 0;

    c[4].vertices[0] = 8;
    c[4].vertices[1] = 5;
    c[4].vertices[2] = 7;
    c[4].vertices[3] = 10;
    c[4].material_id = 0;

    c[5].vertices[0] = 7;
    c[5].vertices[1] = 6;
    c[5].vertices[2] = 9;
    c[5].vertices[3] = 10;
    c[5].material_id = 0;

    c[6].vertices[0] = 8;
    c[6].vertices[1] = 10;
    c[6].vertices[2] = 13;
    c[6].vertices[3] = 11;
    c[6].material_id = 0;

    c[7].vertices[0] = 10;
    c[7].vertices[1] = 9;
    c[7].vertices[2] = 12;
    c[7].vertices[3] = 14;
    c[7].material_id = 0;

    c[8].vertices[0] = 11;
    c[8].vertices[1] = 13;
    c[8].vertices[2] = 15;
    c[8].vertices[3] = 17;
    c[8].material_id = 1;

    c[9].vertices[0] = 13;
    c[9].vertices[1] = 10;
    c[9].vertices[2] = 19;
    c[9].vertices[3] = 15;
    c[9].material_id = 0;

    c[10].vertices[0] = 10;
    c[10].vertices[1] = 14;
    c[10].vertices[2] = 16;
    c[10].vertices[3] = 19;
    c[10].material_id = 0;

    c[11].vertices[0] = 14;
    c[11].vertices[1] = 12;
    c[11].vertices[2] = 18;
    c[11].vertices[3] = 16;
    c[11].material_id = 1;

    c[12].vertices[0] = 17;
    c[12].vertices[1] = 15;
    c[12].vertices[2] = 19;
    c[12].vertices[3] = 20;
    c[12].material_id = 0;

//    c[13].vertices[0] = 15;
//    c[13].vertices[1] = 17;
//    c[13].vertices[2] = 31; //// 31
//    c[13].vertices[3] = 19;
//    c[13].material_id = 1;

    c[13].vertices[0] = 16;
    c[13].vertices[1] = 18;
    c[13].vertices[2] = 21;
    c[13].vertices[3] = 19;
    c[13].material_id = 0;

    c[14].vertices[0] = 20;
    c[14].vertices[1] = 19;
    c[14].vertices[2] = 22;
    c[14].vertices[3] = 23;
    c[14].material_id = 0;

    c[15].vertices[0] = 19;
    c[15].vertices[1] = 21;
    c[15].vertices[2] = 24;
    c[15].vertices[3] = 22;
    c[15].material_id = 0;

    c[16].vertices[0] = 20;
    c[16].vertices[1] = 23;
    c[16].vertices[2] = 26;
    c[16].vertices[3] = 25;
    c[16].material_id = 0;

    c[17].vertices[0] = 23;
    c[17].vertices[1] = 22;
    c[17].vertices[2] = 27;
    c[17].vertices[3] = 26;
    c[17].material_id = 1;

    c[18].vertices[0] = 22;
    c[18].vertices[1] = 24;
    c[18].vertices[2] = 28;
    c[18].vertices[3] = 27;
    c[18].material_id = 1;

    c[19].vertices[0] = 24;
    c[19].vertices[1] = 21;
    c[19].vertices[2] = 29;
    c[19].vertices[3] = 28;
    c[19].material_id = 0;

    printf("%d %d %d %d\n",
    c[0].vertices[0],
    c[0].vertices[1],
    c[0].vertices[2],
    c[0].vertices[3]
    );
    dealii::GridReordering<2>::reorder_cells(c);
    printf("%d %d %d %d\n",
    c[0].vertices[0],
    c[0].vertices[1],
    c[0].vertices[2],
    c[0].vertices[3]
    );
    triangulation .create_triangulation_compatibility (v, c, dealii::SubCellData());
//    triangulation .refine_global (n_ref);

    std::ofstream out ("grid-2.eps");
    dealii::GridOut grid_out;
    grid_out.write_eps (triangulation , out);
};


template <uint8_t dim>
void set_quadrate (dealii::Triangulation<dim> &triangulation, 
        cdbl x0, cdbl x1, cdbl x2, cdbl x3,
        cdbl y0, cdbl y1, cdbl y2, cdbl y3,
        // const double lower, const double top,
        size_t n_refine)
{
    // const double x0 = 0.0;
    // const double x1 = lower;
    // const double x2 = top;
    // const double x3 = 128.0;

//    std::vector< dealii::Point< 2 > > v (8);
//
//    v[0][0] = x0; v[0][1] = x0;
//    v[1][0] = x4; v[1][1] = x0;
//    v[2][0] = x1; v[2][1] = x1;
//    v[3][0] = x3; v[3][1] = x1;
//    v[4][0] = x0; v[4][1] = x4;
//    v[5][0] = x1; v[5][1] = x3;
//    v[6][0] = x3; v[6][1] = x3;
//    v[7][0] = x4; v[7][1] = x4;
//
//    std::vector< dealii::CellData< 2 > > c (5, dealii::CellData<2>());
//
//    c[0].vertices[0] = 0;
//    c[0].vertices[1] = 1;
//    c[0].vertices[2] = 2;
//    c[0].vertices[3] = 3;
//    c[0].material_id = 0;
//
//    c[1].vertices[0] = 0;
//    c[1].vertices[1] = 2;
//    c[1].vertices[2] = 4;
//    c[1].vertices[3] = 5;
//    c[1].material_id = 0;
//    
//    c[2].vertices[0] = 2;
//    c[2].vertices[1] = 3;
//    c[2].vertices[2] = 5;
//    c[2].vertices[3] = 6;
//    c[2].material_id = 1;
//
//    c[3].vertices[0] = 4;
//    c[3].vertices[1] = 5;
//    c[3].vertices[2] = 7;
//    c[3].vertices[3] = 6;
//    c[3].material_id = 0;
//    
//    c[4].vertices[0] = 1;
//    c[4].vertices[1] = 7;
//    c[4].vertices[2] = 3;
//    c[4].vertices[3] = 6;
//    c[4].material_id = 0;

    std::vector< dealii::Point< 2 > > v (16);

    v[0]  = dealii::Point<dim>(x0, y0);
    v[1]  = dealii::Point<dim>(x1, y0);
    v[2]  = dealii::Point<dim>(x2, y0);
    v[3]  = dealii::Point<dim>(x3, y0);
    v[4]  = dealii::Point<dim>(x0, y1);
    v[5]  = dealii::Point<dim>(x1, y1);
    v[6]  = dealii::Point<dim>(x2, y1);
    v[7]  = dealii::Point<dim>(x3, y1);
    v[8]  = dealii::Point<dim>(x0, y2);
    v[9]  = dealii::Point<dim>(x1, y2);
    v[10] = dealii::Point<dim>(x2, y2);
    v[11] = dealii::Point<dim>(x3, y2);
    v[12] = dealii::Point<dim>(x0, y3);
    v[13] = dealii::Point<dim>(x1, y3);
    v[14] = dealii::Point<dim>(x2, y3);
    v[15] = dealii::Point<dim>(x3, y3);

    std::vector< dealii::CellData< 2 > > c (9, dealii::CellData<2>());

    c[6].vertices[0] = 8;  c[7].vertices[0] = 9;  c[8].vertices[0] = 10;
    c[6].vertices[1] = 9;  c[7].vertices[1] = 10; c[8].vertices[1] = 11;
    c[6].vertices[2] = 12; c[7].vertices[2] = 13; c[8].vertices[2] = 14;
    c[6].vertices[3] = 13; c[7].vertices[3] = 14; c[8].vertices[3] = 15;
    c[6].material_id = 0;  c[7].material_id = 0;  c[8].material_id = 0;

    c[3].vertices[0] = 4;  c[4].vertices[0] = 5;  c[5].vertices[0] = 6;
    c[3].vertices[1] = 5;  c[4].vertices[1] = 6;  c[5].vertices[1] = 7;
    c[3].vertices[2] = 8;  c[4].vertices[2] = 9;  c[5].vertices[2] = 10;
    c[3].vertices[3] = 9;  c[4].vertices[3] = 10; c[5].vertices[3] = 11;
    c[3].material_id = 0;  c[4].material_id = 1;  c[5].material_id = 0;

    c[0].vertices[0] = 0;  c[1].vertices[0] = 1;  c[2].vertices[0] = 2;
    c[0].vertices[1] = 1;  c[1].vertices[1] = 2;  c[2].vertices[1] = 3;
    c[0].vertices[2] = 4;  c[1].vertices[2] = 5;  c[2].vertices[2] = 6;
    c[0].vertices[3] = 5;  c[1].vertices[3] = 6;  c[2].vertices[3] = 7;
    c[0].material_id = 0;  c[1].material_id = 0;  c[2].material_id = 0;


    triangulation .create_triangulation (v, c, dealii::SubCellData());

    triangulation .refine_global (n_refine);
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
ATools::FourthOrderTensor unphysical_to_physicaly (
        ATools::FourthOrderTensor &unphys)
{
    enum {x, y, z};
    ATools::FourthOrderTensor res;

    double A = 
        unphys[x][x][x][x] * unphys[y][y][y][y] * unphys[z][z][z][z] - 
        unphys[y][y][z][z] * unphys[z][z][y][y] * unphys[x][x][x][x] +
        unphys[x][x][y][y] * unphys[y][y][z][z] * unphys[z][z][x][x] - 
        unphys[y][y][x][x] * unphys[x][x][y][y] * unphys[z][z][z][z] - 
        unphys[y][y][y][y] * unphys[x][x][z][z] * unphys[z][z][x][x] +
        unphys[y][y][x][x] * unphys[x][x][z][z] * unphys[z][z][y][y]; 

//    printf("%f %f %f A = %f\n", 
//            unphys[x][x][x][x][0], 
//            unphys[y][y][y][y][0], 
//            unphys[z][z][z][z][0], 
//            A);

    for (uint8_t i = 0; i < 3; ++i)
    {
        int no_1 = (i + 1) % 3;
        int no_2 = (i + 2) % 3;

        for (uint8_t j = 0; j < 3; ++j)
        {
            int k = (j == no_1) ? no_2 : no_1;

            if (i == j)
                res[i][i][j][j] = A;
            else
                res[i][i][j][j] = 
                    (unphys[i][i][j][j] * unphys[k][k][k][k] -
                     unphys[i][i][k][k] * unphys[j][j][k][k]);

            res[i][i][j][j] /= 
                (unphys[no_1][no_1][no_1][no_1] * 
                 unphys[no_2][no_2][no_2][no_2] - 
                 unphys[no_1][no_1][no_2][no_2] * 
                 unphys[no_2][no_2][no_1][no_1]);
        };
    };
        
    return res;

};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
void set_circ(dealii::Triangulation< 2 > &triangulation, 
        const double radius, const size_t n_refine)
{
    dealii::GridGenerator ::hyper_cube (triangulation, 0.0, 1.0);
   puts("1111111111111111111111111111111111111111111");
    triangulation .refine_global (n_refine);
    {
        dealii::Point<2> center (0.5, 0.5);
        dealii::Triangulation<2>::active_cell_iterator
            cell = triangulation .begin_active(),
                 end_cell = triangulation .end();
        for (; cell != end_cell; ++cell)
        {
            dealii::Point<2> midle_p(0.0, 0.0);

            for (size_t i = 0; i < 4; ++i)
            {
                midle_p(0) += cell->vertex(i)(0);
                midle_p(1) += cell->vertex(i)(1);
            };
            midle_p(0) /= 4.0;
            midle_p(1) /= 4.0;

           printf("%f %f\n", midle_p(0), midle_p(1));

            if (center.distance(midle_p) < radius)
            {
                cell->set_material_id(1);
//                puts("adf");
            }
            else
                cell->set_material_id(0);
        };
    };
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
void set_cylinder(dealii::Triangulation< 3 > &triangulation, 
        const double radius, cst ort, const size_t n_refine)
{
    dealii::GridGenerator ::hyper_cube (triangulation, 0.0, 1.0);
   // puts("1111111111111111111111111111111111111111111");
    triangulation .refine_global (n_refine);
    {
        dealii::Point<2> center (0.5, 0.5);
        dealii::Triangulation<3>::active_cell_iterator
            cell = triangulation .begin_active(),
                 end_cell = triangulation .end();
        for (; cell != end_cell; ++cell)
        {
            dealii::Point<2> midle_p(0.0, 0.0);

            for (size_t i = 0; i < 8; ++i)
            {
                st count = 0;
                for (st j = 0; j < 3; ++j)
                {
                    if (j != ort)
                    {
                        midle_p(count) += cell->vertex(i)(j);
                        ++count;
                    };
                };
                // midle_p(0) += cell->vertex(i)(0);
                // midle_p(1) += cell->vertex(i)(1);
                // midle_p(0) += cell->vertex(i)(0);
                // midle_p(1) += cell->vertex(i)(2);
            };
            midle_p(0) /= 8.0;
            midle_p(1) /= 8.0;

           // printf("%f %f\n", midle_p(0), midle_p(1));

            if (center.distance(midle_p) < radius)
            {
                cell->set_material_id(1);
//                puts("adf");
            }
            else
                cell->set_material_id(0);
        };
    };
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
void set_ball(dealii::Triangulation< 3 > &triangulation, 
        const double radius, const size_t n_refine)
{
    dealii::GridGenerator ::hyper_cube (triangulation, 0.0, 1.0);
   puts("1111111111111111111111111111111111111111111");
    triangulation .refine_global (n_refine);
    {
        dealii::Point<3> center (0.5, 0.5, 0.5);
        dealii::Triangulation<3>::active_cell_iterator
            cell = triangulation .begin_active(),
                 end_cell = triangulation .end();
        for (; cell != end_cell; ++cell)
        {
            dealii::Point<3> midle_p(0.0, 0.0, 0.0);

            for (size_t i = 0; i < 8; ++i)
            {
                midle_p(0) += cell->vertex(i)(0);
                midle_p(1) += cell->vertex(i)(1);
                midle_p(2) += cell->vertex(i)(2);
            };
            midle_p(0) /= 8.0;
            midle_p(1) /= 8.0;
            midle_p(2) /= 8.0;

           // printf("%f %f\n", midle_p(0), midle_p(1));

            if (center.distance(midle_p) < radius)
            {
                cell->set_material_id(1);
//                puts("adf");
            }
            else
                cell->set_material_id(0);
        };
    };
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
void set_ball(dealii::Triangulation< 3 > &triangulation, 
        const dealii::Point<3> center, cdbl radius)
{
    dealii::Triangulation<3>::active_cell_iterator
        cell = triangulation .begin_active(),
             end_cell = triangulation .end();
    for (; cell != end_cell; ++cell)
    {
        dealii::Point<3> midle_p(0.0, 0.0, 0.0);

        for (size_t i = 0; i < 8; ++i)
        {
            midle_p(0) += cell->vertex(i)(0);
            midle_p(1) += cell->vertex(i)(1);
            midle_p(2) += cell->vertex(i)(2);
        };
        midle_p(0) /= 8.0;
        midle_p(1) /= 8.0;
        midle_p(2) /= 8.0;

        // printf("%f %f\n", midle_p(0), midle_p(1));

        if (center.distance(midle_p) < radius)
        {
            cell->set_material_id(1);
                           // puts("adf");
        }
        // else
        //     cell->set_material_id(0);
        // printf("first %ld\n", cell->material_id());
    };
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
void set_cube(dealii::Triangulation< 3 > &triangulation, 
        const dealii::Point<3> center, cdbl half_len)
{
    dealii::Triangulation<3>::active_cell_iterator
        cell = triangulation .begin_active(),
             end_cell = triangulation .end();
    for (; cell != end_cell; ++cell)
    {
        dealii::Point<3> midle_p(0.0, 0.0, 0.0);

        for (size_t i = 0; i < 8; ++i)
        {
            midle_p(0) += cell->vertex(i)(0);
            midle_p(1) += cell->vertex(i)(1);
            midle_p(2) += cell->vertex(i)(2);
        };
        midle_p(0) /= 8.0;
        midle_p(1) /= 8.0;
        midle_p(2) /= 8.0;

        // printf("%f %f\n", midle_p(0), midle_p(1));

        if (
                (midle_p(0) < (center(0) + half_len)) and 
                (midle_p(1) < (center(1) + half_len)) and 
                (midle_p(2) < (center(2) + half_len)) and 
                (midle_p(0) > (center(0) - half_len)) and 
                (midle_p(1) > (center(1) - half_len)) and 
                (midle_p(2) > (center(2) - half_len))) 
        {
            cell->set_material_id(1);
                           // puts("adf");
        }
        // else
        //     cell->set_material_id(0);
        // printf("first %ld\n", cell->material_id());
    };
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
void set_circ_in_hex(dealii::Triangulation< 2 > &triangulation, 
        const double radius, const size_t n_refine)
{
    cdbl hight = sqrt(3.0);
    dealii::Point<2> p1(0.0, 0.0);
    dealii::Point<2> p2(1.0, hight);
    dealii::GridGenerator ::hyper_rectangle (triangulation, p1, p2);
    triangulation .refine_global (n_refine);
    {
        dealii::Triangulation<2>::active_cell_iterator
            cell = triangulation .begin_active(),
                 end_cell = triangulation .end();
        for (; cell != end_cell; ++cell)
        {
            dealii::Point<2> midle_p(0.0, 0.0);

            for (size_t i = 0; i < 4; ++i)
            {
                midle_p(0) += cell->vertex(i)(0);
                midle_p(1) += cell->vertex(i)(1);
            };
            midle_p(0) /= 4.0;
            midle_p(1) /= 4.0;

            cell->set_material_id(0);
            {
                dealii::Point<2> center (0.5, 0.0);
                if (center.distance(midle_p) < radius)
                    cell->set_material_id(1);
            };
            {
                dealii::Point<2> center (0.0, hight / 2.0);
                if (center.distance(midle_p) < radius)
                    cell->set_material_id(1);
            };
            {
                dealii::Point<2> center (1.0, hight / 2.0);
                if (center.distance(midle_p) < radius)
                    cell->set_material_id(1);
            };
            {
                dealii::Point<2> center (0.5, hight);
                if (center.distance(midle_p) < radius)
                    cell->set_material_id(1);
            };
        };
    };
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
#define SUM(I, BEGIN, END, BODY) ({dbl tmp = 0.0; for (st I = BEGIN; I < END; ++I) {tmp += BODY;}; tmp;});




// template <typename Func, typename Func2>
st foo(cst i, lmbd<st(cst)> &&func, lmbd<st(cst)> &&func2)
{
    return func2(func(i));
};




//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
void set_long_rod (dealii::Triangulation<3> &triangulation, cdbl len, cdbl size, cst n_refine)
{


    dealii::Point<3> p1(0.0, 0.0, 0.0);
    dealii::Point<3> p2(len, 1.0, 1.0);
    dealii::GridGenerator::hyper_rectangle(triangulation, p1, p2);
    triangulation.begin_active()->face(0)->set_boundary_indicator(1);
    triangulation.begin_active()->face(1)->set_boundary_indicator(2);
    triangulation.begin_active()->face(2)->set_boundary_indicator(0);
    triangulation.begin_active()->face(3)->set_boundary_indicator(0);
    triangulation.begin_active()->face(4)->set_boundary_indicator(0);
    triangulation.begin_active()->face(5)->set_boundary_indicator(0);
    triangulation .refine_global(n_refine);
    for (st i = 0; i < 1; ++i)
    {
        dealii::Point<3> center(i * 1.0 + 0.5, 0.5, 0.5);
        // set_cube (triangulation, center, size);
        set_ball (triangulation, center, size);
    };
};

//НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ НЕ ИСПОЛЬЗУЕТСЯ
void set_speciment (dealii::Triangulation<3> &triangulation, 
        cdbl size_x, cdbl size_y, cdbl size_z, cdbl size_inclusion, cdbl size_cell,
        arr<st, 6> id_border, cst n_refine)
{

    dealii::Point<3> p1(0.0, 0.0, 0.0);
    dealii::Point<3> p2(size_x, size_y, size_z);
    dealii::GridGenerator::hyper_rectangle(triangulation, p1, p2);
    triangulation.begin_active()->face(0)->set_boundary_indicator(id_border[0]);
    triangulation.begin_active()->face(1)->set_boundary_indicator(id_border[1]);
    triangulation.begin_active()->face(2)->set_boundary_indicator(id_border[2]);
    triangulation.begin_active()->face(3)->set_boundary_indicator(id_border[3]);
    triangulation.begin_active()->face(4)->set_boundary_indicator(id_border[4]);
    triangulation.begin_active()->face(5)->set_boundary_indicator(id_border[5]);
    triangulation .refine_global(n_refine);
    st num_cell_x = st(size_x / size_cell);
    st num_cell_y = st(size_y / size_cell);
    st num_cell_z = st(size_z / size_cell);
    printf("%ld, %ld %ld\n", num_cell_x, num_cell_y, num_cell_z);
    for (st i = 0; i < num_cell_x; ++i)
    {
        for (st j = 0; j < num_cell_y; ++j)
        {
            for (st k = 0; k < num_cell_z; ++k)
            {
                dealii::Point<3> center((i + 0.5) * size_cell, (j + 0.5) * size_cell, (k + 0.5) * size_cell);
                // set_cube (triangulation, center, size);
                set_ball (triangulation, center, size_inclusion);
            };
        };
    };
};





dbl Analytic_function (const dealii::Point<2> p, cst n, dbl nu)	//посмотреть b(ширину)
{
    cdbl PI = 3.14159265359;
    dbl Uz = 0.0;
    dbl Uw = 0.0;
    dbl b = 2.0;
    dbl c0 = 0.5;
    dbl C0 = 0.125;
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
//    lmbd<st(cst)> add_i = [](cst i){return i-1;};
//    printf("%ld\n", foo(10, [](cst i){return i-1;}, [](cst i){return i+3;}));
//    {
//        Domain<3> domain;
//        {
//            dealii::GridGenerator::hyper_cube(domain.grid, 0.0, 2.0);
//            domain.grid.refine_global(1);
//        };
//    printf("level %d\n", domain.grid.n_global_levels());
//    };

    
    //HEAT_CONDUCTION_NIKOLA_PROBLEM
    if (1)
    {
        Domain<2> domain;
        // {
        //     vec<prmt::Point<2>> boundary_of_segments;
        //     vec<st> types_boundary_segments;
        //     arr<st, 4> types_boundary = {0, 1, 2, 3}; //clockwise
        //     cst num_segments = 1;
        //     prmt::Point<2> p1(0.0, 0.0);
        //     prmt::Point<2> p2(1.0, 1.0);
        //     debputs();
        //     GTools::give_rectangle_with_border_condition (
        //             boundary_of_segments, types_boundary_segments, 
        //             types_boundary, num_segments, p1, p2);
        //     debputs();
        //     make_grid (domain.grid, boundary_of_segments, types_boundary_segments);
        //     domain.grid.refine_global(3);
        //     // for (st i = 0; i < types_boundary_segments.size(); ++i)
        //     // {
        //     //     printf("%ld\n",types_boundary_segments[i]);
        //     // };
        //     // puts("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        // };
		dealii::GridIn<2> gridin;
		gridin.attach_triangulation(domain.grid);
		static char *MeshFileName = "1x2_T1.2.msh";
//		static char *MeshFileName = "1x2_T1.2_nestr.msh";

		std::ifstream f(MeshFileName);


		Assert (dim==2, ExcInternalError());
		gridin.read_msh(f);
	    domain.grid.refine_global(5);
        	std::cout << "\tNumber of active cells:       "
                	<< domain.grid.n_active_cells()
                	<< std::endl;

        debputs();
        dealii::FE_Q<2> fe(1);
        domain.dof_init (fe);

        SystemsLinearAlgebraicEquations slae;
        ATools ::trivial_prepare_system_equations (slae, domain);


        cdbl c0 = 0.5;
		cdbl C0 = 0.125;
        cdbl E = 1.0;
        cdbl nu = 0.25;
        cdbl mu = E / ( 2 * ( 1 + nu ) );
//		cdbl mu = 1.0;
        LaplacianScalar<2> element_matrix (domain.dof_handler.get_fe());
        {
//            element_matrix.C .resize(1);
//            element_matrix.C[0][x][x] = mu;
//            element_matrix.C[0][x][y] = 0.0;
//            element_matrix.C[0][y][x] = 0.0;
//            element_matrix.C[0][y][y] = mu;

            element_matrix.C .resize(2);
            element_matrix.C[0][x][x] = mu;
            element_matrix.C[0][x][y] = 0.0;
            element_matrix.C[0][y][x] = 0.0;
            element_matrix.C[0][y][y] = mu;
            element_matrix.C[1][x][x] = mu;
            element_matrix.C[1][x][y] = 0.0;
            element_matrix.C[1][y][x] = 0.0;
            element_matrix.C[1][y][y] = mu;
        };

//        vec<arr<typename Nikola::SourceScalar<2>::Func, 2>> U(1);
//        U[0][x] = [mu, nu, c0] (const dealii::Point<2> &p) {return mu*nu*0.5*(  (p(0)-c0)*(p(0)-c0) - p(1)*p(1)  );}; //Ux
//        U[0][y] = [mu, nu, c0] (const dealii::Point<2> &p) {return mu*nu*(p(0)-c0)*p(1);}; //Uy
        vec<arr<typename Nikola::SourceScalar<2>::Func, 2>> U(2);
        U[0][x] = [mu, nu, c0, C0] (const dealii::Point<2> &p) {return mu*nu*0.5*(std::pow(p(0)-c0,2.0)-std::pow(p(1),2.0)+2.0*C0);}; //Ux
        U[0][y] = [mu, nu, c0] (const dealii::Point<2> &p) {return mu*nu*(p(0)-c0)*p(1);}; //Uy
        U[1][x] = [mu, nu, c0, C0] (const dealii::Point<2> &p) {return mu*nu*0.5*(std::pow(p(0)-c0,2.0)-std::pow(p(1),2.0)+2.0*C0);}; //Ux
        U[1][y] = [mu, nu, c0] (const dealii::Point<2> &p) {return mu*nu*(p(0)-c0)*p(1);};


//        vec<typename Nikola::SourceScalar<2>::Func> tau(1);
//        tau[0] = [E, c0] (const dealii::Point<2> &p) {return E*(p(0)-c0);};		//умножить на -1 и записать здесь

        vec<typename Nikola::SourceScalar<2>::Func> tau(2);
        tau[0] = [E, c0] (const dealii::Point<2> &p) {return E*(p(0)-c0);};
        tau[1] = [E, c0] (const dealii::Point<2> &p) {return E*(p(0)-c0);};


        // vec<arr<typename Nikola::SourceScalar<2>::Func, 2>> U(2);
        // U[0][x] = [] (const dealii::Point<2> &p) {return 0.0;}; //Ux
        // U[0][y] = [] (const dealii::Point<2> &p) {return 0.0;}; //Uy
        // U[1][x] = [] (const dealii::Point<2> &p) {return 0.0;};
        // U[1][y] = [] (const dealii::Point<2> &p) {return 0.0;};
        // vec<typename Nikola::SourceScalar<2>::Func> tau(2);
        // tau[0] = [] (const dealii::Point<2> &p) {return -2.0;};
        // tau[1] = [] (const dealii::Point<2> &p) {return -2.0;};
        Nikola::SourceScalar<2> element_rhsv (U, tau, domain.dof_handler.get_fe());

        Assembler::assemble_matrix<2> (slae.matrix, element_matrix, domain.dof_handler);
        Assembler::assemble_rhsv<2> (slae.rhsv, element_rhsv, domain.dof_handler);

        


        HCPTools ::print_temperature<2> (slae.rhsv, domain.dof_handler, "b");
        dbl sum = 0.0;
        for (st i = 0; i < slae.rhsv.size(); ++i)
        {
            sum += slae.rhsv(i);
        };
        printf("Integral %f\n", sum);

//	vec<BoundaryValueScalar<2>> bound (1);
//        bound[0].function      = [] (const dealii::Point<2> &p) {return 0;}; // функция граничного условия
//	bound[0].boundary_id   = 0; // номер границ к которым применяется это условие
//	bound[0].boundary_type = TBV::Neumann; // Тип граничного условия, в данном случае Дирихле
//	for (auto b : bound)
//	    ATools ::apply_boundary_value_scalar<2> (b) .to_slae (slae, domain);
////	    std::cout << "make boundary\n";


        dealii::SolverControl solver_control (10000, 1e-12);
        dealii::SolverCG<> solver (solver_control);
        solver.solve (
                slae.matrix,
                slae.solution,
                slae.rhsv
                ,dealii::PreconditionIdentity()
                );

//=============================================================================
    FILE *FileOfResults;
    FileOfResults = fopen ("FileOfResults.txt","a");
//					«r»		-	Режим открытия файла для чтения. Файл должен существовать.
//					«w»		-	Режим создания пустого файла для записи. Если файл с таким
//								именем уже существует его содержимое стирается, и файл
//								рассматривается как новый пустой файл.
//					«a»		-	Дописать в файл. Операция добавления данных в конец файла.
//								Файл создается, если он не существует.
//					«r+»	-	Режим открытия фала для обновления чтения и записи.
//								Этот файл должен существовать.
//					«w+»	-	Создаёт пустой файл для чтения и записи. Если файл с таким
//								именем уже существует его содержимое стирается, и файл
//								рассматривается как новый пустой файл.
//					«a+»	-	Открыть файл для чтения и добавления данных. Все операции
//								записи выполняются в конец файла, защищая предыдущее
//								содержания файла от случайного изменения. Вы можете
//								изменить позицию (FSEEK, перемотка назад) внутреннего
//								указателя на любое место файла только для чтения, операции
//								записи будет перемещать указатель в конец файла, и только
//								после этого дописывать новую информацию. Файл создается,
//								если он не существует.

		dbl LimitOfDifference = 0.05;
		dbl MaxDiff = -1000;
		dbl MaxDiff_otn = -1000;
        dealii::Vector<dbl> analytic(slae.solution.size());
        dealii::Vector<dbl> diff(slae.solution.size());
        dealii::Vector<dbl> diff_otn(slae.solution.size());
        for (auto cell = domain.dof_handler.begin_active (); cell != domain.dof_handler.end (); ++cell)
        {
            for (st i = 0; i < dealii::GeometryInfo<2>::vertices_per_cell; ++i)
            {
                dbl indx = cell->vertex_dof_index(i, 0);
                analytic(indx) = Analytic_function(cell->vertex(i), 32, nu);
				if( std::abs( analytic( indx ) ) >= LimitOfDifference)
					{
			            diff(indx) = std::abs(analytic(indx) - slae.solution(indx));
						if (MaxDiff < diff(indx)) {MaxDiff = diff(indx);}
			            diff_otn(indx) = diff(indx) / (slae.solution(indx) ) * 100.0;
						if (MaxDiff_otn < diff_otn(indx)) {MaxDiff_otn = diff_otn(indx);}
					}
            };
        };

        HCPTools ::print_temperature<2> (slae.solution, domain.dof_handler, "temperature.gpd");
        HCPTools ::print_temperature<2> (analytic, domain.dof_handler, "Analytic_function.gpd");
        HCPTools ::print_temperature<2> (diff, domain.dof_handler, "diff.gpd");
        HCPTools ::print_temperature<2> (diff_otn, domain.dof_handler, "diff_otn.gpd");
        printf("MaxDiff 	= %f\n", MaxDiff);
        printf("MaxDiff_otn	= %f %%\n", MaxDiff_otn);
		
		fprintf(FileOfResults,"MeshFileName 	= %s\n", MeshFileName);
//		fprintf(FileOfResults,"MeshFileName 	= \n");
		fprintf(FileOfResults,"\tNumber of active cells:       %d\n", domain.grid.n_active_cells());
		fprintf(FileOfResults,"MaxDiff 	= %f\n", MaxDiff);
		fprintf(FileOfResults,"MaxDiff_otn 	= %f %%\n\n\n\n", MaxDiff_otn);

		fclose(FileOfResults);	
    };

    return EXIT_SUCCESS;

}


//    FILE *FileName;
//    FileName = fopen ("FileName.txt","w");
//					«r»		-	Режим открытия файла для чтения. Файл должен существовать.
//					«w»		-	Режим создания пустого файла для записи. Если файл с таким
//								именем уже существует его содержимое стирается, и файл
//								рассматривается как новый пустой файл.
//					«a»		-	Дописать в файл. Операция добавления данных в конец файла.
//								Файл создается, если он не существует.
//					«r+»	-	Режим открытия фала для обновления чтения и записи.
//								Этот файл должен существовать.
//					«w+»	-	Создаёт пустой файл для чтения и записи. Если файл с таким
//								именем уже существует его содержимое стирается, и файл
//								рассматривается как новый пустой файл.
//					«a+»	-	Открыть файл для чтения и добавления данных. Все операции
//								записи выполняются в конец файла, защищая предыдущее
//								содержания файла от случайного изменения. Вы можете
//								изменить позицию (FSEEK, перемотка назад) внутреннего
//								указателя на любое место файла только для чтения, операции
//								записи будет перемещать указатель в конец файла, и только
//								после этого дописывать новую информацию. Файл создается,
//								если он не существует.
//
//    fprintf(FileName,"  Fa\t\t\t   F\t\t\t  abserror\t\t    mape\n");
//    fprintf(FileName,"\t\tMaxAbs:  %8.4e\t\tMaxMape:  %8.4e\n", MaxAbs, MaxMape);
//
//
//    fclose(FileName);

