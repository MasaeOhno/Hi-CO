# Hi-CO
3D genome structure analysis with nucleosome resolution 

## Contents
- LabVIEW code for analyzing paired-end sequencing reads to output the Hi-CO matrix
- Sample paired-end sequencing data (R1 and R2) with SAM format
- FORTRAN code for creating an initial structure of SA-MD simulation

## Installation and usage

Hi-CO matrix 
1. Install LabVIEW and Vision Development Module from National Instruments.
2. Download "Create_Hi-COmatrix".
3. Input names of the Read1 and Read2 files (.sam) to "SAM file path 1" and "SAM file path 2".
4. Input a list of nucleosome positions to "nucleosome positions (chr)" and "nucleosome positions (base)", and select "Edit/Make current values default"
5. Input a list of nucleosome numbers in every chromosome to "nucleosome number", and a list of base pair length of every chromosome to "chromosome length base", and select "Edit/Make current values default".
6. Run.

  
Initial structure for simulation
- Install gfortran.
- Download "Generate_initial_model".
- The following command line is used to generate an array of the nucleosome particles. `<N>` is the number of nucleosomes.
```
gfortran genini.f90 -o genini.exe
```
```
./genini.exe <N>
```
- The putput is "out.pqr", which can be visualized by VMD, UCSF chimera, etc.
