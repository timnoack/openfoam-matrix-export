#include <fstream>

#include "error.H"
#include "fast_matrix_market/fast_matrix_market.hpp"
#include "fast_matrix_market/types.hpp"
#include "matrixExporter.H"

namespace Foam {
defineTypeNameAndDebug(matrixExporter, 0);

lduMatrix::solver::addsymMatrixConstructorToTable<matrixExporter>
    addsmoothSolverSymMatrixConstructorToTable_;

lduMatrix::solver::addasymMatrixConstructorToTable<matrixExporter>
    addsmoothSolverAsymMatrixConstructorToTable_;
}  // namespace Foam

Foam::matrixExporter::matrixExporter(
    const word &fieldName, const lduMatrix &matrix,
    const FieldField<Field, scalar> &interfaceBouCoeffs,
    const FieldField<Field, scalar> &interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList &interfaces,
    const dictionary &solverControls)
    : lduMatrix::solver(fieldName, matrix, interfaceBouCoeffs,
                        interfaceIntCoeffs, interfaces, solverControls) {
  readControls();
}

void Foam::matrixExporter::readControls() {
  lduMatrix::solver::readControls();
  fileNameA_ = controlDict_.getOrDefault<fileName>("matrixFile", "A.mtx");
  fileNameb_ = controlDict_.getOrDefault<fileName>("sourceFile", "b.mtx");
  fileNamex_ = controlDict_.getOrDefault<fileName>("initialFieldFile", "x.mtx");
  comment_ =
      controlDict_.getOrDefault<string>("comment", "No description provided");
  exitAfterExport_ = controlDict_.getOrDefault<bool>("exitAfterExport", true);
}

void Foam::matrixExporter::exportMatrix() const {
  Info << "Exporting matrix to file " << fileNameA_ << endl;

  label numRows = matrix_.lduAddr().size();
  label numFaces = matrix_.lduAddr().lowerAddr().size();
  label numNonZeros = numRows + 2 * numFaces;

  fast_matrix_market::matrix_market_header header;
  header.comment = comment_;
  header.nrows = numRows;
  header.ncols = numRows;
  header.format = fast_matrix_market::coordinate;
  header.field = fast_matrix_market::real;
  header.symmetry = matrix_.symmetric() ? fast_matrix_market::symmetric
                                        : fast_matrix_market::general;

  std::vector<label> rows;
  std::vector<label> cols;
  std::vector<scalar> values;

  rows.reserve(numNonZeros);
  cols.reserve(numNonZeros);
  values.reserve(numNonZeros);

  auto addCoeffs = [&](const label row, const label col, const scalar coeff) {
    rows.push_back(row);
    cols.push_back(col);
    values.push_back(coeff);
  };

  forAll(matrix_.lduAddr(), row) {
    // Add diagonal coefficient
    addCoeffs(row, row, matrix_.diag()[row]);

    // Add lower coefficients of this row
    label ownerStart = matrix_.lduAddr().ownerStartAddr()[row];
    label ownerEnd = matrix_.lduAddr().ownerStartAddr()[row + 1];
    for (label facei = ownerStart; facei < ownerEnd; ++facei) {
      label neighbourCell = matrix_.lduAddr().upperAddr()[facei];
      scalar coeff = matrix_.lower()[facei];
      addCoeffs(row, neighbourCell, coeff);
    }

    // Add upper coefficients of this row if matrix is not symmetric
    // If the matrix is symmetric, the matrix is marked as symmetric in the
    // matrix market header and the upper coefficients MUST NOT be added
    if (matrix_.symmetric()) continue;

    label losortStart = matrix_.lduAddr().losortStartAddr()[row];
    label losortEnd = matrix_.lduAddr().losortStartAddr()[row + 1];
    for (label i = losortStart; i < losortEnd; ++i) {
      label facei = matrix_.lduAddr().losortAddr()[i];
      label ownerCell = matrix_.lduAddr().lowerAddr()[facei];
      scalar coeff = matrix_.upper()[facei];
      addCoeffs(row, ownerCell, coeff);
    }
  }

  std::ofstream os(fileNameA_.c_str());

  if (!os) {
    FatalErrorIn("matrixExporter::exportMatrix()")
        << "Cannot open file " << fileNameA_ << endl;
  }

  fast_matrix_market::write_matrix_market_triplet(os, header, rows, cols,
                                                  values);
  os.close();
}

void Foam::matrixExporter::exportField(const scalarField &field,
                                       fileName fileName) const {
  Info << "Exporting field to file " << fileName << endl;
  label numRows = matrix_.lduAddr().size();

  fast_matrix_market::matrix_market_header header;
  header.comment = comment_;
  header.nrows = numRows;
  header.ncols = 1;
  header.format = fast_matrix_market::array;
  header.field = fast_matrix_market::real;

  std::vector<scalar> values;
  values.reserve(numRows);

  forAll(field, row) { values.push_back(field[row]); }

  std::ofstream os(fileName.c_str());

  if (!os) {
    FatalErrorIn("matrixExporter::field()")
        << "Cannot open file " << fileName << endl;
  }

  fast_matrix_market::write_matrix_market_array(os, header, values);
  os.close();
}

Foam::solverPerformance Foam::matrixExporter::solve(
    scalarField &psi_s, const scalarField &source, const direction cmpt) const {
  solverPerformance solverPerf(typeName, fieldName_);

  Info << "Export-solver invoked for field " << fieldName_ << endl;

  if (UPstream::parRun()) {
    FatalErrorIn("matrixExporter::solve()")
        << "Parallel run not supported. Run the solver application with mpirun."
        << endl;
  }

  if (!interfaces_.empty()) {
    // Interfaces should only exist in parallel runs, but lets make sure anyway
    FatalErrorIn("matrixExporter::solve()")
        << "Interfaces not supported" << endl;
  }

  // Export matrix, source and initial field
  exportMatrix();
  exportField(source, fileNameb_);
  exportField(psi_s, fileNamex_);

  Info << "Export-solver finished for field " << fieldName_ << endl;

  if (exitAfterExport_) {
    error e(
        "Exiting after successfull matrix export. To disable this, set "
        "exitAfterExport to no in the solver settings.");
    e.exit(0);
  }

  return solverPerf;
}
