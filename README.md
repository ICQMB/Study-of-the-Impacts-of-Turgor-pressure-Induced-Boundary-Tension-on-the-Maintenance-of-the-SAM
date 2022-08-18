# ScePlantCells
This repository contains C++ code for modeling plant cells using the subcellular element model.

The code was compiled in the linux environment 4.18.0-348.12.2.el8_5.x86_64
using the g++ (GCC)  8.5.0 20210514 (Red Hat 8.5.0-4) compiler.

Compile by moving all .cpp and .h files in the root directory of this repository into the same directory and execute
"make clean" // Removes pre-existing build and object files
"make all" // Builds new executable.

Executables can be run via ./program -option1 -option2...

The first six options MUST be ordered as follows:
-option1 - path to output folder for VTK animations with cytoplasm rendered
-option2 - path to output folder for node-scale output data
-option3 - path to output folder for nematic output data
-option4 - path to output folder for VTK animations without cytoplasm rendered
-option5 - path to output folder for cell-scale data
-option6 - path to output folder for tisue-scale data

All subsequent options have the form
-key [value]

All keys are found in main.cpp, but the main key-value pairs used for this paper are...

-OOP_flag [1 or 0] // 1 is P3D model; 0 is 2D model
-BFM [double X] // sets X micronewtons as the magnitude of F_boundary.

This repository includes several script files used to expedite the submission process in a SLURM queue system.
