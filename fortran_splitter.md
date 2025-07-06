# Fortran Code Splitter

A Python utility to split legacy Fortran code into multiple files based on a mapping configuration.

## Features

- Splits a single Fortran source file into multiple files
- Preserves original formatting and comments
- Automatically handles dependencies (COMMON blocks, PARAMETERs)
- Warns about subroutine/function call dependencies
- Separates main program into its own file (`*_main.f90`)
- Case-insensitive matching for Fortran identifiers

## Usage

```bash
python fortran_splitter.py source_file.f90 mapping_file.txt
```

## Mapping File Format

The mapping file uses a simple YAML-like format:

```
output_file1.f90:
- SUBROUTINE1
- FUNCTION1

output_file2.f90:
- SUBROUTINE2
- FUNCTION2
```

## Example

Given a source file `gridv.f90` and mapping file `split_config.txt`:

```bash
python fortran_splitter.py gridv.f90 split_config.txt
```

This will generate:
- `gridv_main.f90` (containing the main PROGRAM)
- `gridv_io.f90` (containing I/O-related routines)
- `gridv_calculations.f90` (containing calculation routines)
- etc.

## How It Works

1. **Parsing**: The utility parses both the mapping file and Fortran source
2. **Analysis**: It identifies:
   - Main program
   - Subroutines and functions
   - COMMON blocks and PARAMETER statements
   - Dependencies between routines
3. **Generation**: For each output file, it:
   - Includes only the COMMON blocks and PARAMETERs used by routines in that file
   - Adds the specified subroutines/functions
   - Preserves original formatting

## Warnings

The utility will warn you about:
- Routines that call other routines in different files
- Routines specified in the mapping file but not found in the source
- Unmapped routine dependencies

## Requirements

- Python 3.6+
- No external dependencies (uses only standard library)

## Limitations

- Supports Fortran 77 and Fortran 90/95 style code
- Does not handle modules (beyond USE statements)
- Does not reorder routines to resolve dependencies