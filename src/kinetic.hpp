#pragma once
#ifndef DFT_KINETIC_HPP
#define DFT_KINETIC_HPP
#include "fespace.h"
#include <iostream>

class KineticOperator
{
  public:
    explicit KineticOperator(const DFTGLLHexSpace &space);

    void Assemble();

    bool IsAssembled() const;

    const mfem::SparseMatrix &Matrix() const;

    void Mult(const mfem::Vector &x_true, mfem::Vector &y_true) const;

    double Energy(const mfem::Vector &x_true) const;

    void PrintInfo(std::ostream &os = std::cout) const;

  private:
    static mfem::IntegrationRule make_tensor_gll_rule_(int order);

    const DFTGLLHexSpace &m_space;
    mfem::FiniteElementSpace *m_fes;
    int m_order;

    mfem::SparseMatrix m_matrix;
    bool m_is_assembled{false};

    mfem::ConstantCoefficient m_half_coefficient{0.5};
};
#endif // DFT_KINETIC_HPP