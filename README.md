# MPScode

This is an experimental program to test and apply DMRG methods on medium scale systems with complex entangelement structures.
It uses variatonal MPS-based DMRG for the computation of ground states and low-lying excited states of one-dimensional quantum
systems. The utilized algorithm is the DMRG3S [1], supporting the use abelian symmetries to greatly reduce the computational
effort. Next to ground- and excited state search, an iDMRG implementation is contained, too, which targets ground states
of infinite systems based on two-site DMRG. Additional support for two-site DMRG exists, but is currently not interfaced
outside the iDMRG.

The program internally uses arpack++ [2] to iteratively solve the eigenvalue problem which has to be solved in the 
matrix optimization step.

[1] C.Hubig et al. Phys. Rev. B 91, 155115
[2] https://www.ime.unicamp.br/~chico/arpack++/
