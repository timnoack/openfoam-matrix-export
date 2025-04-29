#ifndef Foam_matrixExport_exportSolver_H
#define Foam_matrixExport_exportSolver_H
#include "lduMatrix.H"
#include "runTimeSelectionTables.H"

namespace Foam {

class matrixExporter : public lduMatrix::solver {
private:
    fileName destDir_;
    string comment_;
    bool exitAfterExport_;
    autoPtr<lduMatrix::solver> actualSolver_;

protected:
    //- Read the control parameters from the controlDict_
    virtual void readControls();

    fileName getPathTo(word type) const;

    void exportMatrix() const;
    void exportField(const scalarField &field, word type) const;

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

} // namespace Foam
#endif