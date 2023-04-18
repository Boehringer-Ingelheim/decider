## R CMD check results

0 errors | 0 warnings | 3 notes

* This is a new release.

* Note: GNU make is a SystemRequirements. 
  Explanation: This is due to the use of the rstan package for MCMC sampling
  and Bayesian computations.
  
* Note containing the following information:
  \code{
       Found 'abort', possibly from 'abort' (C), 'runtime' (Fortran)
       Found 'exit', possibly from 'exit' (C), 'stop' (Fortran)
       Found 'printf', possibly from 'printf' (C)
    
    Compiled code should not call entry points which might terminate R nor
    write to stdout/stderr instead of to the console, nor use Fortran I/O
    nor system RNGs. The detected symbols are linked into the code but
    might come from libraries and not actually be called.}
    
  Explanation: The package contains C++ code that was generated using the rstan 
  package, which seems to contain the relevant "abort"/"exit"/"printf" 
  statements triggering the above note. None of these statements are called in
  the functions included in this package or by the included stanmodels.

