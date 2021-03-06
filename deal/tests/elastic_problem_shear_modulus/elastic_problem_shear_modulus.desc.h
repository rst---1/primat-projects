/*
 * =====================================================================================
 *
 *       Filename:  heat_conduction_problem.desc.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  14.09.2012 10:55:04
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef HEAT_CONDUCTION_PROBLEM_DESC

#define HEAT_CONDUCTION_PROBLEM_DESC

#include <array>
#include <projects/deal/tests/esm_elastic_problem/esm_elastic_problem.h>
#include <projects/deal/main/problem/problem.h>
#include <projects/deal/main/solver_se/solver_se.h>
#include <projects/deal/main/function/function.h>
#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/matrices.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/solver_cg.h>
#include <fstream>
#include <iostream>


template <int dim>
class BoundaryValues : public dealii::Function<dim>
{
  public:
    BoundaryValues ();

    virtual void vector_value (const dealii::Point<dim> &p,
                               dealii::Vector<double>   &values) const;

    void set (typename ElasticProblemSup<dim>::TypeFunc &value);

    typename ElasticProblemSup<dim>::TypeFunc content;
};

template<uint8_t dim>
class ElementRightHandSideVectorElasticProblemShearModulus : 
    public ElementRightHandSideVector<
    dim, double, typename ElasticProblemSup<dim>::TypeFunc>
{
    public:
        ElementRightHandSideVectorElasticProblemShearModulus ();

    virtual void set_coefficient (typename ElasticProblemSup<dim>::TypeFunc &coef); 

    virtual double operator() (const size_t index_i, 
            const dealii::QGauss<dim> &quadrature_formula, 
            const dealii::FEValues<dim> &fe_values) const;
};

template <uint8_t dim>
class ElasticProblemShearModulus : public Problem< 
                              dim, 
                              ElementStiffnessMatrixElasticProblem<dim>, 
                              ElementRightHandSideVectorElasticProblem<dim> >
{
    public:
        ElasticProblemShearModulus (const dealii::Triangulation<dim> &triangulation,
                        typename ElasticProblemSup<dim>::TypeCoef &coefficient, 
                        typename ElasticProblemSup<dim>::TypeFunc &boundary_values,
                        typename ElasticProblemSup<dim>::TypeFunc &rhs_values);
        ~ElasticProblemShearModulus ();
    //Methods
    public:
        virtual Report solved ();
        virtual void print_result (const std::string &filename);

    protected:
        virtual Report setup_system ();
        virtual Report assemble_system ();
        virtual Report apply_boundary_values ();
        virtual Report solve_system_equations ();
        virtual Report calculate_shear_modulus ();
        virtual Report output_results ();

    //Fields
    protected:
        BoundaryValues<dim>   boundary_values;
        dealii::FESystem<dim> finite_element;
        std::string           output_file_name;
        double area_of_domain;
        double meta_shear_modulus;
};

#endif
