# tdse_basis
Solve Time-dependent Schrodinger Equation (TDSE) using numerical basis in Hilbert Space. 

The program is based on MPI and PETSc libraries. It runs on multiple nodes of a computer cluster.

npsflib-basis: This directory includes source code of core functions. It is built as a library and will be used by a main program.

example: This directory includes basic examples of main program.

ais: This directory includes an example to calculate attosecond transient absoprtion for autoinoization states of He atom. 

ne_absorption: This directory includes an example to calculate attosecond transient absoprtion for Ne atom. 

Reference: Phys. Rev. A 86, 013410 (2012).
