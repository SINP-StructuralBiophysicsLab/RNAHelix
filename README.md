[![License](https://img.shields.io/badge/License-GPL3.0-blue.svg?style=flat-square)](https://github.com/Manisso/fsociety/blob/master/LICENSE) ![OS](https://img.shields.io/badge/Tested%20On-Linux%20|%20Windows-yellowgreen.svg?style=flat-square) 

# RNAHelix: version 1.0
```bash

  ____     _       _                _     _   _______    _          _  _       _  
 |  _ \   | \     | |     /\       | |   | |  |  ___ |  | |        | | \ \    / / 
 | | \ \  | |\    | |     \ \      | |   | |  | |   \|  | |        | |  \ \  / /  
 | |  | | | | \   | |    / \ \     | |___| |  | |__     | |        | |   \ \/ /   
 | |_/ /  | |  \  | |   / / \ \    | |___| |  |  __|    | |        | |    \ \/    
 |  _ /   | |   \ | |  / /___\ \   | |   | |  | |       | |        | |   / \ \   
 | | \ \  | |    \| | / /_____\ \  | |   | |  | |___/|  | |____/|  | |  / / \ \  
 |_|  \_\ |_|     \_|/ /       \_\ |_|   |_|  |______|  |_______|  |_| /_/   \_\  
===================================================================================
```
                                -----=====%-+-+-+%=====-----

## || Development and Citation ||
### RNAHelix:

- Bhattacharyya,D., Halder,S., Basu,S., Mukherjee,D., Kumar,P. and Bansal,M. (2017)
RNAHelix: computational modeling of nucleic acid structures with Watson–Crick and 
noncanonical base pairs. J. Comput. Aided. Mol. Des., 31, 219–235.
### NUCGEN:

- S. Mukherjee, M. Bansal and D. Bhattacharyya (2006) Conformational 
   specificity of non-canonical base pairs 4 and higher order structures 
   in nucleic acids: crystal 5 structure database analysis, J. Comp. 
   Aided Mol. Des. 20, 629-45. 
   
   
- Bansal,M., Bhattacharyya,D. and Ravi,B. (1995) NUPARM and NUCGEN: software 
   for analysis and generation of sequence dependent nucleic acid structures.
   Comput. Appl. Biosci., 11, 281-287.

                                -----=====%-+-+-+%=====-----

## || Contents ||

       |>RNAHelix
	 |
	 --DataSet.dat
	 |
	 --prm2param.pl
	 |
         --parameter.loc
         |
         --RNAHelix.f
         |
	 --RNAHelix.exe
	 |
	 --RNAHelix.linux
	 |
	 --par_all27_na.inp
	 |
	 --top_all27_rna.inp
	 |
	 --top_all27_dna.inp
	 |
	 --minimize.inp
	 |
	 --README
	 |
	 --LICENSE
 	 |
	 --sample/


If you do not receive all the files mentiond above then there is some problem
  Download the.tar again to run it properly.
  
                                -----=====%-+-+-+%=====-----

### ||Requirements||

 G77 compiler or Intel ifort compiler is required for running the programme.
 
                                -----=====%-+-+-+%=====-----

# Installation

## Installation [Linux](https://wikipedia.org/wiki/Linux) [![alt tag](http://icons.iconarchive.com/icons/dakirby309/simply-styled/32/OS-Linux-icon.png)](https://fr.wikipedia.org/wiki/Linux)



 - Global:
  ------
```bash
   The precompiled version requires the file "DataSet.dat" to be placed in
   /usr/local/bin . In case you are unable to copy the supplied DataSet.dat 
   file into /usr/local/bin, please modify the source-code (line no. 687) to
   specify the correct location of the file and recompile using a suitable
   FORTRAN-77 compiler.
```
```bash
  -----
   For local installation just put the uncompressed file in desired location 
   and run the software from any where using the full path or you can add the
   folder location to your PATH variable.
   if you have put your folder in /home/Username/bin/RnaHelix, for example,	
   then
    bash:
        PATH=$PATH:/home/Username/bin/RnaHelix
	export PATH
    tcsh:
	setenv PATH $PATH\:/home/Username/bin/RnaHelix
```
                 -----=====%-+-+-+%=====-----

## ||Run programme||
```bash
  ./RNAHelix.g77 


Input [mandatory]			Description
-------------------------------------------------------------------------------
  DataSet.dat                   ideal parameter files for nucleotides with 
				Watson-Crick, Sugar and Hoogsteeien edges this
				file is important for the generation of nucleotides

  parameter.loc			rotational, translational and vibrational 
				parameters to be used for the generation of 
				nucleic acid structure. This file is prepared 
                                by processing the nuparm output *.prm with prm2param.pl


Output [default]			Description
-------------------------------------------------------------------------------
  allatoms.pdb			the regenerated nucleotide structure file in .pdb
				format

  centers.bpc		        basepair centers

  double.hlx                    basepairing information of the regenarated helix
```

                                -----=====%-+-+-+%=====-----


## || File Descriptions ||

```bash
  
DataSetSugar.dat
<><><><><><><><>
  
  this file is segmented into 15 segments (four nucleotides[ACGT/U] with all 
  three edge of binding Watson-Crick, Sugar, and Hoosteeine). It contains the 
  ideal parameters for helix generation and, therefore, SHOULD NOT BE MODIFIED.

parameter.loc
<><><><><><><>
  
  1st line
  ---------------
  col 1-4 	int	number of base-pair	
  col 8		0/1	0 for no sequence file and 1 for using sequence file
  col 12
  col 16-		name of the prm file used for generation
  subsequent lines
  -----------------
  col 1-3		nucleotide information
  col 5-8		basepairing type
  col 13-18		Tilt
  col 21-26		Roll
  col 29-34		Twist
  col 37-42		Shift
  col 45-50		Slide
  col 53-58		Rise
  col 61-66		Buckle
  col 69-74       	Open
  col 77-82       	Propel
  col 85-90       	Stagger
  col 93-98       	Shear
  col 101-106     	Stretch
  
                Tilt    Roll   Twist    Shift  Slide   Rise    Buckle  Open    Propel   Stagger Shear Stretch
  C:G W:WC     -0.34    5.69   29.18    0.20   -2.09    3.30    6.15   -1.10   -5.85   -0.03    0.31    2.94
  
  
allatoms.pdb
<><><><><><>

  this file is the primary output of the RnaHelixGen software. this file is is 
  in format of pdb file. one can visualize the file in Pymol/Rasmol/VMD or any 
  standard visualizing tools.

centers.bpc
<><><><><>
 
double.hlx
<><><><><>

  this file shows a textual visualization of the structure with cannonical and/or 
  non-cannonical basepairing information, where '|' stands for a Watson-Crick edge
  pairing, '-' stands for Sugar edge pairing and '*' stands for Hoosteeine edge 
  pairing

Additional files for Backbone Generation by CHARMM
<><><><><><>
   Please use one of the supplied topology files appropriate for your molecule.
   Minimize the allatoms.pdb file using minimize.inp script, which reads the
   sequence of the DNA/RNA from Sequence.str file.
<><><><><><>

We recommend use of BPFIND and NUPARM as additional data preparation tools.
These also can be downloaded from [Bioinformatics](http://www.saha.ac.in/biop/bioinformatics.html) page.
```
