# OpenFOAM matrix exporter
This library registers the `matrixExporter` solver in OpenFOAM. When the solver is invoked, it will export the matrix it is supposed to solve to a file in the Matrix Market format.

## Usage
1. Clone the library to any folder. You may clone it into the `$WM_PROJECT_USER_DIR` folder, but you can use any other folder as well.
```bash
git clone git@github.com:timnoack/openfoam-matrix-export.git
```
2. Build the library using `wmake libso`:
```bash
wmake libso path/to/openfoam-matrix-export
```
3. Add the library to the `system/controlDict` of your case:
```bash
libs ("libMatrixExporter");
```
4. Use the `matrixExporter` solver in your `system/fvSolution` file:
```bash
solvers
{
    T // the name of your field
    {
        solver           matrixExporter;
        matrixFile       "A.mtx";
        sourceFile       "b.mtx";
        initialFieldFile "x.mtx";
        exitAfterExport  yes;
        comment          "Please provide a short description of the testcase and its matrix here.";
    }
}
```
5. Run your case as usual. The solver will export the matrix, source vector and initial field to the specified filenames.