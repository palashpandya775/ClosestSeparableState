# ClosestSeparableState
Closest separable state to a given reference state using the Gilbert's Algorithm for optimizing a quadratic form over a convex set.
(Current code only looks for the closest fully separable state. To be updated later.)

The file 'gilbert.c' is the main code file. The program can be compiled using the makefile, or using the following command:

    gcc gilbert.c -o a.out -lgsl -lgslcblas -lm -O3

*Reguires GNU Scientific Library*

Input: Reference state in a file with the naming convention 'ppt<#>.dat'
Output: Two files- 
        1. Containing the approximations to the closest separable state: "pptm<#>*d*-*n*r.txt"
        2. Containing the sequence of Hilbert-Schmidt distances:  "pptm<#>*d*-*n*.txt"

The command to run the program takes 4 integers:
1. Number of parties, n
2. Dimension of the subsystem, d
3. Number in the input file, <#>
   (Also appears in the names of the output files.)
4. Number of Corrections, cs

  
An example:
  
The provided input file 'ppt1.dat', contains the 2 qubit Bell state in the format that is read using gsl_matrix_fscanf.
To run the program for this state we compile the program and run the output file:
  
      $-> a.out 2 2  1   1000  
                n d <#>   cs  
  
The output is in the corresponding files: pptm12-2.txt and pptm12-2r.txt.
    
Another Example:
    
The file 'ppt2.dat', contains the 3 qubit GHZ state, and the program is run  by:
      
      $-> a.out 3 2  1   1000  
                n d <#>   cs 
  
