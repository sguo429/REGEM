# REGEM

REGEM (GEM-reanalysis-tool) is a software program for re-analyzing large-scale gene-environment interaction testing results, including multi-exposure interaction, joint, and marginal tests. It uses results directly from [GEM](https://github.com/large-scale-gxe-methods/GEM) output.

<br />
Current version: 1.0

<br />

## Contents
- [Quick Installation](#quick-installation)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Contact](#contact)
- [License](#license)

## Quick Installation

To install REGEM, run the following lines of code:
 ```
git clone https://github.com/large-scale-gxe-methods/REGEM
cd REGEM
cd src
make
 ```
<br />
<br />
<br />

## Dependencies

- boost-1.70.0 or up
- Intel Math Kernel Library: tested with mkl-2019.3.199

## Usage

### Running REGEM

1. [Command Line Options](#command-line-options)
1. [Input Files](#input-files)
1. [Output File Format](#output-file-format)
1. [Examples](#examples)

<br /> 
<br />

### Command Line Options
Once REGEM is compiled, the executable ./REGEM can be used to run the program.
For a list of options, use```cd REGEM/src./REGEM -help```.

<details>
     <summary> <b>List of Options</b> </summary>

```
General Options:

   --help                
	Prints available options and exits.
   --version             
	Prints the version of REGEM and exits.


Input File Options:
   --results-file        
	Path to the GEM results file.
   --out                
	Full path and extension to where REGEM output results.
        Default: regem.out


Phenotype File Options:
   --exposure-names      
	One or more column names in the phenotype file naming the exposure(s) to be included in interaction tests.
   --int-covar-names     
	Any column names in the phenotype file naming the covariate(s) for which interactions should be included for adjustment (mutually exclusive with --exposure-names).


Filtering Options:
   --maf                 
	Threshold to filter variants based on the minor allele frequency.
        Default: 0.001
   --miss-geno-cutoff    
	Threshold to filter variants based on the missing genotype rate.
        Default: 0.05
```
</details>

<br /> 

### Input Files

REGEM directly uses output files from GEM (v1.4.1 and up).

### Output File Format

REGEM will write results to the output file specified with the --out parameter (or 'regem.out' if no output file is specified).
Below are details of the possible column headers in the output file:

```diff 
SNPID              - The SNP identifier as retrieved from the genotype file.
CHR                - The chromosome of the SNP.
POS                - The physical position of the SNP. 
Non_Effect_Allele  - The allele not counted in association testing.  
Effect_Allele      - The allele that is counted in association testing.  
N_Samples          - The number of samples without missing genotypes.
AF                 - The allele frequency of the effect allele.      

Beta_Marginal           - The coefficient estimate for the marginal genetic effect (i.e., from a model with no interaction terms).
SE_Beta_Marginal        - The model-based SE associated with the marginal genetic effect estimate.  
robust_Beta_Marginal  
robust_SE_Beta_Marginal - The robust SE associated with the marginal genetic effect estimate.

Beta_G             - The coefficient estimate for the genetic main effect (G).
Beta_G-*           - The coefficient estimate for the interaction or interaction covariate terms.
SE_Beta_G          - Model-based SE associated with the the genetic main effect (G).  
SE_Beta_G-*        - Model-based SE associated with any GxE or interaction covariate terms.
Cov_Beta_G_G-*          - Model-based covariance between the genetic main effect (G) and any GxE or interaction covariate terms.  
robust_Beta_G  
robust_Beta_G-*    
robust_SE_Beta_G   - Robust SE associated with the the genetic main effect (G).  
robust_SE_Beta_G-* - Robust SE associated with any GxE or interaction covariate terms.
robust_Cov_Beta_G_G-*   - Robust covariance between the genetic main effect (G) and any GxE or interaction covariate terms.   

P_Value_Marginal           - Marginal genetic effect p-value from model-based SE.
P_Value_Interaction        - Interaction effect p-value (K degrees of freedom test of interaction effect) from model-based SE. (K is number of major exposures)
P_Value_Joint              - Joint test p-value (K+1 degrees of freedom test of genetic and interaction effect) from model-based SE.
robust_P_Value_Marginal    - Marginal genetic effect p-value from robust SE.
robust_P_Value_Interaction - Interaction effect p-value from robust SE.
robust_P_Value_Joint       - Joint test p-value (K+1 degrees of freedom test of genetic and interaction effect) from robust SE.
```

<br />

### Examples
<br />

```unix
./REGEM --results-file GEM.out --exposure-names cov1 --int-covar-names cov2 --out regem.out
```
<br />
<br />

## Contact 
For comments, suggestions, bug reports and questions, please contact Han Chen (Han.Chen.2@uth.tmc.edu), Alisa Manning (AKMANNING@mgh.harvard.edu), Kenny Westerman (KEWESTERMAN@mgh.harvard.edu) or Cong Pan (Cong.Pan@uth.tmc.edu). For bug reports, please include an example to reproduce the problem without having to access your confidential data.

<br />
<br />

## License 

 ```
 GEM : Gene-Environment interaction analysis for Millions of samples
 Copyright (C) 2018-2021  Liang Hong, Han Chen, Duy Pham, Cong Pan
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ```
 
