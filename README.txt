-----------------------------------------------------------
 How to use the derivative-free optimizer DFL for MINLP
-----------------------------------------------------------

1- Gunzip and untar the archive in a folder on your computer.

2- FORTRAN90 version of the code

   Edit file problem.f90 to setup your own problem.
   In particular, modify the subroutines 
   setdim    : which sets problem dimension
   setbounds : which sets upper and lower bounds on the variables
   startp    : which sets the starting point
   funct     : which defines the objective function

   PYTHON version of the code
   Edit file PYTHON/problem.py to setup your own problem.
   In particular, modify the procedures 
   setbounds : which sets upper and lower bounds on the variables
   startp    : which sets the starting point
   funct     : which defines the objective function

3- At command prompt execute 

     $> make
 
   which will create the executable 'sdbox'

5- execute

     $> ./dfbox

   or

     $> python PYTHON/main_box.py
