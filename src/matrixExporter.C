#include <fstream>

#include "Time.H"
#include "dictionary.H"
#include "error.H"
#include "fast_matrix_market/fast_matrix_market.hpp"
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
                        interfaceIntCoeffs, interfaces, solverControls),
      actualSolver_(nullptr) {
  readControls();
}

void Foam::matrixExporter::readControls() {
  lduMatrix::solver::readControls();
  destDir_ =
      controlDict_.getOrDefault<fileName>("directory", "./matrixExport/");
  comment_ =
      controlDict_.getOrDefault<string>("comment", "No description provided");

  if (dictionary *solverConfig = controlDict_.findDict("solverConfig")) {
    actualSolver_ =
        lduMatrix::solver::New(fieldName_, matrix_, interfaceBouCoeffs_,
                               interfaceIntCoeffs_, interfaces_, *solverConfig);
    Info << "Export-Solver uses actual solver " << actualSolver_->type()
         << endl;
  }

  // Exit after export by default if no actual solver is provided
  exitAfterExport_ = controlDict_.getOrDefault<bool>("exitAfterExport",
                                                     actualSolver_ == nullptr);
}

Foam::fileName Foam::matrixExporter::getPathTo(word type) const {
  const Time &time = matrix_.mesh().thisDb().time();
  word timeName = time.timeName();

  fileName currentDir = destDir_ / timeName;

  if (!exists(currentDir)) mkDir(currentDir);

  return currentDir / fieldName_ + "_" + type + ".mtx";
}

void Foam::matrixExporter::exportMatrix() const {
  fileName outputFileA = getPathTo("matrix");
  Info << "Exporting matrix to file " << outputFileA << endl;

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

  std::ofstream os(outputFileA.c_str());

  if (!os) {
    FatalErrorIn("matrixExporter::exportMatrix()")
        << "Cannot open file " << outputFileA << endl;
  }

  fast_matrix_market::write_matrix_market_triplet(os, header, rows, cols,
                                                  values);
  os.close();
}

void Foam::matrixExporter::exportField(const scalarField &field,
                                       word type) const {
  fileName outputFileName = getPathTo(type);
  Info << "Exporting field to file " << outputFileName << endl;
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

  std::ofstream os(outputFileName.c_str());

  if (!os) {
    FatalErrorIn("matrixExporter::field()")
        << "Cannot open file " << outputFileName << endl;
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
        << "Parallel run not supported. Run the non-decomposed case on a "
           "single core."
        << endl;
  }

  if (!interfaces_.empty()) {
    // Interfaces should only exist in parallel runs, but lets make sure anyway
    FatalErrorIn("matrixExporter::solve()")
        << "Interfaces not supported" << endl;
  }

  // Export matrix, source and initial field
  exportMatrix();
  exportField(source, "source");
  exportField(psi_s, "psi_initial");

  Info << "Export-solver finished for field " << fieldName_ << endl;

  // If actualSolver is available, use the configured solver to solve the matrix
  if (actualSolver_) {
    Info << "Solving matrix for field " << fieldName_ << " using solver "
         << actualSolver_->type() << endl;

    // Solve the matrix
    solverPerf = actualSolver_->solve(psi_s, source, cmpt);

    // Export the solution field
    exportField(psi_s, "psi_solution.mtx");
  }

  if (exitAfterExport_) {
    error e(
        "Exiting after successful matrix export. To disable this, set "
        "exitAfterExport to no in the solver settings.");
    e.exit(0);
  }

  return solverPerf;
}