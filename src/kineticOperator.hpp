#pragma once

#include "FEspace.h"
#include <iostream>

class KineticOperator
{
  public:
    explicit KineticOperator(const DFTGLLHexSpace &space);

    void Assemble();

    bool IsAssembled() const { return assembled_; }

    const mfem::SparseMatrix &Matrix() const;

    void Mult(const mfem::Vector &x_true, mfem::Vector &y_true) const;

    double Energy(const mfem::Vector &x_true) const;

    void PrintInfo(std::ostream &os = std::cout) const;

  private:
    static mfem::IntegrationRule MakeTensorGLLRule_(int order);

    const DFTGLLHexSpace &space_;
    mfem::FiniteElementSpace *fes_;
    int order_;

    mfem::SparseMatrix matrix_;
    bool assembled_{false};

    mfem::ConstantCoefficient half_coeff_{0.5};
};
