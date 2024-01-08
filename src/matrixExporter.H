#ifndef Foam_matrixExport_exportSolver_H
#define Foam_matrixExport_exportSolver_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

class matrixExporter : public lduMatrix::solver {
 private:
  fileName fileNameA_;
  fileName fileNameb_;
  fileName fileNamex_;
  string comment_;

 protected:
  // Protected Data
  //- Read the control parameters from the controlDict_
  virtual void readControls();

  void exportMatrix() const;
  void exportField(const scalarField &field, fileName fileName) const;

 public:
  //- Runtime type information
  TypeName("matrixExporter");

  // Constructors

  //- Construct from matrix components and solver controls
  matrixExporter(const word &fieldName, const lduMatrix &matrix,
                 const FieldField<Field, scalar> &interfaceBouCoeffs,
                 const FieldField<Field, scalar> &interfaceIntCoeffs,
                 const lduInterfaceFieldPtrsList &interfaces,
                 const dictionary &solverControls);

  //- Destructor
  virtual ~matrixExporter() = default;

  // Member Functions

  //- Solve the matrix with this solver
  virtual solverPerformance solve(scalarField &psi, const scalarField &source,
                                  const direction cmpt = 0) const;
};

}  // namespace Foam

#endif
