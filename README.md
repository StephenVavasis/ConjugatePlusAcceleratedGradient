# ConjugatePlusAcceleratedGradient
## Hybrid conjugate gradient/accelerated gradient for minimizing smooth convex functions

The Matlab function `cplusag` in this repository minimizes a smooth convex function *f* via nonlinear conjugate gradient.  A sufficient-progress measure
defined in terms of estimating sequences from Nesterov, *Introductory Lectures on Convex Optimization*, Kluwer 2004, is invoked on each iteration, and if
insufficient progress is made, the method switches to accelerated gradient.  Refer to our paper Karimi & Vavasis, Nonlinear conjugate gradient for smooth
convex functions, available on arxiv.org, 2021.

A version of accelerated gradient is available in `ag`.  In addition, there are three smooth convex test cases available.  Function `makeorthol1test` computes
the ABPDN (approximate basis-pursuit denoising).  Function `makehandles_logistic` computes logistic regression, and `make_huber_regression`
computes Huber regerssion.  All three test problems are described in the paper.

Hager & Zhang's CG-Descent is available here: http://users.clas.ufl.edu/hager/papers/Software/  We have used version 6.8.  Note that version 6.8 sometimes crashes
in Matlab because it runs out of memory.  This repository contains a modified version of the Matlab front-end, called `cg_descent.c`, that is more aggressive
in freeing up memory from intermediate results and seems to fix the out-of-memory issues.  Note that the file `cg_descent.c` here is meant to
replace the file `cg_descent.c` that appears in Hager & Zhang's subdirectory called MATLAB, not the file with the same name in the main
directory.

