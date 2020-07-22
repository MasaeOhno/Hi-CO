# Hi-CO
3D genome structure analysis with nucleosome resolution 

Contents
- LabVIEW code for analyzing paired-end sequencing reads to output the Hi-CO matrix
- Sample paired-end sequencing data (R1 and R2) with SAM format
- FORTRAN code for creating an initial structure of SA-MD simulation

Installation & usage

-Hi-CO matrix 
  Install LabVIEW and Vision Development Module from National Instruments.
  Download "Create_Hi-COmatrix".
  Input locations and names of the Read1 and Read2 files (.sam) to "HiC sam file path 1" and "HiC sam file path 2".
  Run
  
-Initial structure for simulation
  Install gfortran.
  Download "Generate_Model.zip" and unzip.
  The following command line is used to generate an array of the nucleosome particles.
  1) gfortran genini.f90 -o genini.exe
  2) ./genini.exe <N>
	where <N> is the number of nucleosomes
  The putput is "out.pqr", which can be visualized by VMD, UCSF chimera, etc.
