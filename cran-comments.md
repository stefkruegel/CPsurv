## Resubmission
This is a resubmission. 
In this version I have changed the examples so that there's no conflict with CRAN policies regarding run time


## Test environments
* local Windows x86-64, R 3.2.3
* Debian 3.11.8-1, R 3.1.0


## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

** running examples for arch 'i386' ... [174s] NOTE
Examples with CPU or elapsed time > 10s
             user system elapsed
plot.cpsurv 87.95   0.04   88.05
CPsurv      83.83   0.03   83.87

** running examples for arch 'x64' ... [169s] NOTE
Examples with CPU or elapsed time > 10s
             user system elapsed
plot.cpsurv 84.28   0.02   84.39
CPsurv      81.75   0.06   81.86

 The function can be executed parallelized for shorter run time, but it's disabled in the examples.


## Downstream dependencies
I have also checked downstream dependencies with devtools:revdep_check()
