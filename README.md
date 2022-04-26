# ClosestSeparableState
Closest separable state to a given reference state using the Gilbert's Algorithm for optimizing a quadratic form over a convex set.

The file 'gilbert.c' is the main code file. The program can be compiled using the makefile, or using the following command:

    gcc gilbert.c -o a.out -lgsl -lgslcblas -lm -O3

*Reguires GNU Scientific Library*

Input: Reference state in a file with the naming convention 'ppt<#>.dat'
Output: Two files- 
        1. Containing the approximations to the closest separable state: "pptm<#><d>-<n>r.txt"
        2. Containing the sequence of Hilbert-Schmidt distances:  "pptm<#><d>-<n>.txt"

The command to run the program takes 4 integers:
1. Dimension of the subsystem, d
2. Number of parties, n
3. Number of Corrections, cs
4. Number in the input file, <#>
   (Also appears in the names of the output files.)

  
An example:
  
The provided input file 'ppt1.dat', contains the 2 qubit Bell state in the format that is read using gsl_matrix_fscanf.
To run the program for this state we compile the program and run the output file:
  
      $-> a.out 2 2 1000  1 
                d n  cs  <#> 
  
The output is in the corresponding files: pptm12-2.txt and pptm12-2r.txt.
