#!/usr/bin/env python3
"""
Fortran Code Splitter Utility
Splits a Fortran source file into multiple files based on a mapping file.
"""

import re
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional
import argparse


class FortranSplitter:
    def __init__(self, source_file: str, mapping_file: str):
        self.source_file = source_file
        self.mapping_file = mapping_file
        self.source_content = ""
        self.subroutines = {}  # name -> (start_line, end_line, content)
        self.functions = {}    # name -> (start_line, end_line, content)
        self.common_blocks = []  # List of (line_num, content)
        self.parameters = []     # List of (line_num, content)
        self.module_uses = []    # List of (line_num, content)
        self.main_program = None # (start_line, end_line, content)
        self.file_mapping = {}   # target_file -> [routine_names]
        self.dependencies = defaultdict(set)  # routine -> set of called routines
        self.routine_globals = defaultdict(dict)  # routine -> {'commons': set, 'parameters': set}
        
    def parse_mapping_file(self):
        """Parse the mapping file to get file->subroutines mapping"""
        current_file = None
        
        with open(self.mapping_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                    
                if line.endswith(':'):
                    # New file section
                    current_file = line[:-1].strip()
                    self.file_mapping[current_file] = []
                elif line.startswith('- '):
                    # Subroutine/function name
                    if current_file:
                        routine_name = line[2:].strip()
                        self.file_mapping[current_file].append(routine_name)
                        
    def parse_fortran_source(self):
        """Parse Fortran source to extract components"""
        with open(self.source_file, 'r') as f:
            self.source_content = f.read()
            
        lines = self.source_content.split('\n')
        
        # Extract global declarations first
        self._extract_globals(lines)
        
        # Extract program units
        self._extract_program_units(lines)
        
        # Analyze dependencies
        self._analyze_dependencies()
        
    def _extract_globals(self, lines: List[str]):
        """Extract COMMON blocks, PARAMETERs, and USE statements"""
        in_program_unit = False
        
        for i, line in enumerate(lines):
            line_upper = line.upper().strip()
            
            # Check if we're entering a program unit
            if (re.match(r'^\s*(PROGRAM|SUBROUTINE|FUNCTION|.*FUNCTION)\s+', line_upper) or
                re.match(r'^\s*PROGRAM\s+', line_upper)):
                in_program_unit = True
                continue
                
            # Check if we're exiting a program unit
            if line_upper.startswith('END PROGRAM') or line_upper.startswith('END SUBROUTINE') or \
               line_upper.startswith('END FUNCTION') or line_upper == 'END':
                in_program_unit = False
                continue
                
            # Only collect globals outside of program units
            if not in_program_unit:
                if 'COMMON' in line_upper and not line_upper.startswith('!') and not line_upper.startswith('C'):
                    self.common_blocks.append((i, line))
                elif 'PARAMETER' in line_upper and not line_upper.startswith('!') and not line_upper.startswith('C'):
                    self.parameters.append((i, line))
                elif line_upper.startswith('USE '):
                    self.module_uses.append((i, line))
                    
    def _extract_program_units(self, lines: List[str]):
        """Extract main program, subroutines, and functions"""
        i = 0
        while i < len(lines):
            line = lines[i]
            line_upper = line.upper().strip()
            
            # Skip comments and empty lines
            if not line_upper or line_upper.startswith('!') or line_upper.startswith('C'):
                i += 1
                continue
                
            # Check for PROGRAM
            if re.match(r'^\s*PROGRAM\s+(\w+)', line_upper):
                match = re.match(r'^\s*PROGRAM\s+(\w+)', line_upper)
                prog_name = match.group(1)
                start = i
                end = self._find_end(lines, i, 'PROGRAM')
                self.main_program = (start, end, '\n'.join(lines[start:end+1]))
                i = end + 1
                continue
                
            # Check for SUBROUTINE
            if re.match(r'^\s*SUBROUTINE\s+(\w+)', line_upper):
                match = re.match(r'^\s*SUBROUTINE\s+(\w+)', line_upper)
                sub_name = match.group(1)
                start = i
                end = self._find_end(lines, i, 'SUBROUTINE')
                content = '\n'.join(lines[start:end+1])
                self.subroutines[sub_name] = (start, end, content)
                # Extract local globals for this subroutine
                self._extract_routine_globals(sub_name, lines[start:end+1])
                i = end + 1
                continue
                
            # Check for FUNCTION (including type declarations like DOUBLE PRECISION FUNCTION)
            func_match = re.match(r'^\s*(?:.*\s+)?FUNCTION\s+(\w+)', line_upper)
            if func_match:
                func_name = func_match.group(1)
                start = i
                end = self._find_end(lines, i, 'FUNCTION')
                content = '\n'.join(lines[start:end+1])
                self.functions[func_name] = (start, end, content)
                # Extract local globals for this function
                self._extract_routine_globals(func_name, lines[start:end+1])
                i = end + 1
                continue
                
            i += 1
            
    def _extract_routine_globals(self, routine_name: str, routine_lines: List[str]):
        """Extract COMMON blocks and PARAMETERs used within a routine"""
        self.routine_globals[routine_name] = {'commons': set(), 'parameters': set()}
        
        for line in routine_lines:
            line_upper = line.upper().strip()
            if not line_upper or line_upper.startswith('!') or line_upper.startswith('C'):
                continue
                
            # Extract COMMON block names
            if 'COMMON' in line_upper:
                # Match patterns like COMMON /BLOCKNAME/ or COMMON /BLOCKNAME/VAR1,VAR2
                common_matches = re.findall(r'/(\w+)/', line)
                for block in common_matches:
                    self.routine_globals[routine_name]['commons'].add(block.upper())
                    
            # Check if routine uses PARAMETER variables
            if 'PARAMETER' in line_upper:
                self.routine_globals[routine_name]['parameters'].add('USES_PARAMETERS')
                
    def _find_end(self, lines: List[str], start: int, unit_type: str) -> int:
        """Find the END statement for a program unit"""
        i = start + 1
        while i < len(lines):
            line_upper = lines[i].upper().strip()
            if line_upper.startswith(f'END {unit_type}') or \
               (line_upper == 'END' and self._is_end_of_unit(lines, i, start)):
                return i
            i += 1
        return len(lines) - 1
        
    def _is_end_of_unit(self, lines: List[str], end_idx: int, start_idx: int) -> bool:
        """Check if a bare END statement belongs to the current unit"""
        # Simple heuristic: count nested structures
        depth = 0
        for i in range(start_idx + 1, end_idx):
            line_upper = lines[i].upper().strip()
            if re.match(r'^\s*(SUBROUTINE|FUNCTION|PROGRAM)\s+', line_upper):
                depth += 1
            elif re.match(r'^\s*END\s+(SUBROUTINE|FUNCTION|PROGRAM)', line_upper):
                depth -= 1
        return depth == 0
        
    def _analyze_dependencies(self):
        """Analyze which routines call which other routines"""
        all_routines = set(self.subroutines.keys()) | set(self.functions.keys())
        
        # Check each routine for calls to other routines
        for routine_name, (start, end, content) in {**self.subroutines, **self.functions}.items():
            for other_routine in all_routines:
                if other_routine != routine_name:
                    # Case-insensitive search for CALL statements or function calls
                    call_pattern = rf'\bCALL\s+{other_routine}\b'
                    func_pattern = rf'\b{other_routine}\s*\('
                    
                    if re.search(call_pattern, content, re.IGNORECASE) or \
                       re.search(func_pattern, content, re.IGNORECASE):
                        self.dependencies[routine_name].add(other_routine)
                        
    def generate_output_files(self):
        """Generate the output files based on the mapping"""
        source_path = Path(self.source_file)
        source_stem = source_path.stem
        
        # Collect warnings
        warnings = []
        
        # Generate main program file if it exists
        if self.main_program:
            main_file = f"{source_stem}_main.f90"
            self._write_main_file(main_file)
            print(f"Generated: {main_file}")
            
        # Generate mapped files
        for target_file, routine_names in self.file_mapping.items():
            content_parts = []
            
            # Add file header comment
            content_parts.append(f"! Generated from {self.source_file}")
            content_parts.append(f"! Contains: {', '.join(routine_names)}")
            content_parts.append("")
            
            # Collect needed globals
            needed_commons = set()
            needed_parameters = False
            
            for routine_name in routine_names:
                routine_name_upper = routine_name.upper()
                if routine_name_upper in self.routine_globals:
                    needed_commons.update(self.routine_globals[routine_name_upper]['commons'])
                    if 'USES_PARAMETERS' in self.routine_globals[routine_name_upper]['parameters']:
                        needed_parameters = True
                        
            # Add module uses if any
            if self.module_uses:
                for _, use_stmt in self.module_uses:
                    content_parts.append(use_stmt)
                content_parts.append("")
                
            # Add needed COMMON blocks
            for _, common_line in self.common_blocks:
                # Check if this common block is needed
                common_blocks_in_line = set(re.findall(r'/(\w+)/', common_line.upper()))
                if common_blocks_in_line & needed_commons:
                    content_parts.append(common_line)
                    
            if needed_commons:
                content_parts.append("")
                
            # Add PARAMETER statements if needed
            if needed_parameters:
                for _, param_line in self.parameters:
                    content_parts.append(param_line)
                content_parts.append("")
                
            # Add routines
            for routine_name in routine_names:
                routine_name_upper = routine_name.upper()
                
                # Check dependencies
                if routine_name_upper in self.dependencies:
                    deps = self.dependencies[routine_name_upper]
                    # Check if any dependencies are not in the same file
                    missing_deps = []
                    for dep in deps:
                        dep_found = False
                        for fname, routines in self.file_mapping.items():
                            if dep in [r.upper() for r in routines]:
                                if fname != target_file:
                                    missing_deps.append((dep, fname))
                                dep_found = True
                                break
                        if not dep_found:
                            missing_deps.append((dep, "NOT MAPPED"))
                            
                    if missing_deps:
                        for dep, location in missing_deps:
                            warnings.append(f"WARNING: {routine_name} calls {dep} which is in {location}")
                
                # Add the routine
                if routine_name_upper in self.subroutines:
                    _, _, content = self.subroutines[routine_name_upper]
                    content_parts.append(content)
                elif routine_name_upper in self.functions:
                    _, _, content = self.functions[routine_name_upper]
                    content_parts.append(content)
                else:
                    warnings.append(f"WARNING: Routine {routine_name} not found in source file")
                    
                content_parts.append("")  # Empty line between routines
                
            # Write the file
            with open(target_file, 'w') as f:
                f.write('\n'.join(content_parts))
                
            print(f"Generated: {target_file}")
            
        # Print warnings
        if warnings:
            print("\n" + "="*60)
            print("WARNINGS:")
            print("="*60)
            for warning in warnings:
                print(warning)
                
    def _write_main_file(self, filename: str):
        """Write the main program file"""
        content_parts = []
        
        # Add file header
        content_parts.append(f"! Generated from {self.source_file}")
        content_parts.append("! Main program")
        content_parts.append("")
        
        # Add module uses
        if self.module_uses:
            for _, use_stmt in self.module_uses:
                content_parts.append(use_stmt)
            content_parts.append("")
            
        # Add COMMON blocks
        if self.common_blocks:
            for _, common_line in self.common_blocks:
                content_parts.append(common_line)
            content_parts.append("")
            
        # Add PARAMETERs
        if self.parameters:
            for _, param_line in self.parameters:
                content_parts.append(param_line)
            content_parts.append("")
            
        # Add main program
        _, _, main_content = self.main_program
        content_parts.append(main_content)
        
        with open(filename, 'w') as f:
            f.write('\n'.join(content_parts))


def main():
    parser = argparse.ArgumentParser(description='Split Fortran source file into multiple files')
    parser.add_argument('source_file', help='Fortran source file to split')
    parser.add_argument('mapping_file', help='Mapping file specifying how to split the source')
    
    args = parser.parse_args()
    
    # Validate input files exist
    if not Path(args.source_file).exists():
        print(f"Error: Source file '{args.source_file}' not found")
        sys.exit(1)
        
    if not Path(args.mapping_file).exists():
        print(f"Error: Mapping file '{args.mapping_file}' not found")
        sys.exit(1)
        
    # Create splitter and process files
    splitter = FortranSplitter(args.source_file, args.mapping_file)
    
    print(f"Parsing mapping file: {args.mapping_file}")
    splitter.parse_mapping_file()
    
    print(f"Parsing Fortran source: {args.source_file}")
    splitter.parse_fortran_source()
    
    print("\nGenerating output files...")
    splitter.generate_output_files()
    
    print("\nDone!")


if __name__ == "__main__":
    main()