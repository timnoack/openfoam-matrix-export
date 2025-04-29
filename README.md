# OpenFOAM Matrix and Field Exporter

This library registers the `matrixExporter` solver in OpenFOAM. When the solver is invoked, it exports both the system matrix and the associated fields (source and initial solution field) to files in the Matrix Market format. If configured with a solver, it can also solve the system and export the final solution field. All exports are organized by simulation time and field name.

## Usage

1. Clone the library to any folder. You may clone it into the `$WM_PROJECT_USER_DIR` folder, but you can use any other folder as well. Don't forget to clone the submodules as well (`--recursive`) as the library depends on the [fast_matrix_market](https://github.com/alugowski/fast_matrix_market.git).
```bash
git clone --recursive https://github.com/timnoack/openfoam-matrix-export.git
```

2. Build the library. The shared library is automatically copied into the `$FOAM_USER_LIBBIN` folder.
```bash
wmake libso path/to/openfoam-matrix-export
```

3. Load the library to your solver application by specifying it in your `system/controlDict`:
```cpp
libs ("libMatrixExporter");
```

4. Use the `matrixExporter` solver in your `system/fvSolution` file:
```cpp
solvers
{
    p // the name of your field
    {
        solver          matrixExporter;
        directory       "./matrixExport/"; // Base directory for exported files
        comment         "Please provide a short description of the testcase and its matrix here."; // Optional comment
        
        // Optional but highly recommended: Add a solver configuration to actually solve the matrix
        solverConfig
        {
            solver          PCG;
            preconditioner  DIC;
            tolerance       1e-6;
            relTol          0.01;
            maxIter         1000;
        }
    }
}
```

5. Run your case as usual. 

## Features

### Matrix and Field Exports
Every time the solver is invoked, it exports:
- The complete system matrix (A)
- The source field (b in Ax=b)
- The initial field values (initial x in Ax=b)
- The solution field (if a solver is configured)

All files are written in the Matrix Market format:
- System matrices use the coordinate format for sparse storage
- Fields use the array format for dense vector storage

### File Organization
The matrix and field files are organized in a systematic directory structure:
```
<directory>/<time>/<fieldName>_<type>.mtx
```
Where:
- `<directory>` is the base directory specified in the configuration
- `<time>` is the current simulation time (e.g., "0.5", "1", etc.)
- `<fieldName>` is the name of the field being solved (e.g., "p", "U", "T")
- `<type>` identifies the content:
  - `matrix` - The system matrix (A in Ax=b)
  - `source` - The right-hand side vector (b in Ax=b)
  - `psi_initial` - The initial field values (initial x)
  - `psi_solution` - The solved field values (only if a solver is configured)

### Solver Integration
The exporter can optionally solve the matrix after exporting it:

1. Add a `solverConfig` dictionary to your solver setup
2. The exporter will use that solver configuration to solve the matrix
3. Both the initial field and the final solution will be exported
4. By default, the simulation will continue after exporting and solving the matrix. You can change this behavior by explicitly setting `exitAfterExport` to `yes`.

### Exit Behavior
- By default, if no solver is configured, the exporter will exit after export (`exitAfterExport = yes`)
- If a solver is configured, it will default to continuing the simulation (`exitAfterExport = no`)
- You can override this behavior by explicitly setting `exitAfterExport`

⚠️ **CRITICAL WARNING**: If you set `exitAfterExport = no` without providing a `solverConfig`, the simulation will continue with an unsolved matrix. This means:

- The field values will remain at their initial values
- These incorrect values will be used in calculations for other fields
- Physics coupling between fields will be broken
- In iterative algorithms (SIMPLE, PIMPLE), convergence will be impossible
- Final results will be physically meaningless
- The simulation may appear to run normally, but ALL RESULTS WILL BE INCORRECT

## Examples

### Basic Export (No Solving)

For exporting matrices and fields without solving:

```cpp
solvers
{
    p
    {
        solver          matrixExporter;
        directory       "./matrixExport/";
        comment         "Pressure equation from cavity case";
        // exitAfterExport defaults to yes when no solver is configured
    }
}
```

This will produce these files for each timestep:
```
./matrixExport/0.1/p_matrix.mtx
./matrixExport/0.1/p_source.mtx
./matrixExport/0.1/p_psi_initial.mtx
```

### Export with Solving

For exporting matrices and fields, then solving and exporting the solution:

```cpp
solvers
{
    p
    {
        solver          matrixExporter;
        directory       "./matrixExport/";
        comment         "Pressure equation from cavity case";
        // exitAfterExport defaults to no when a solver is configured
        
        solverConfig
        {
            solver          PCG;
            preconditioner  DIC;
            tolerance       1e-6;
            relTol          0.01;
            maxIter         1000;
        }
    }
}
```

This will produce these files for each timestep:
```
./matrixExport/0.1/p_matrix.mtx  (The system matrix)
./matrixExport/0.1/p_source.mtx  (The right-hand side)
./matrixExport/0.1/p_psi_initial.mtx  (The initial field values)
./matrixExport/0.1/p_psi_solution.mtx  (The solved field values)
```

The workflow:
1. Export the system matrix, source vector, and initial field
2. Solve the pressure equation using PCG with DIC preconditioner
3. Export the solution field
4. Continue the simulation with the solved values

## Working with Exported Files

The exported files are in the Matrix Market format which can be read by many scientific computing tools:

- **MATLAB/Octave**: `A = mmread('p_matrix.mtx'); b = mmread('p_source.mtx');`
- **Python with SciPy**: `from scipy.io import mmread; A = mmread('p_matrix.mtx')`
- **Julia**: `using MatrixMarket; A = MatrixMarket.mmread('p_matrix.mtx')`

This allows for external analysis, visualization, or solving of OpenFOAM's linear systems.