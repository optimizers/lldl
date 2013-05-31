# LLDL

## Overview

LLDL is a modification of [Lin and Moré's limited-memory Cholesky factorization](http://dx.doi.org/10.1137/S1064827597327334) code-named ICFS for symmetric positive definite matrices. LLDL implements a similar limited-memory scheme for symmetric indefinite matrices that possess a **LDL**<sup>T</sup> factorization, i.e., with **D** diagonal. Symmetric quasi-definite matrices fall into this category.

The main idea is that if **L** and **D** are the exact factors of **A**, the preconditioned matrix (**L** |**D**| **L**<sup>T</sup>)<sup>-1</sup> **A** has only two eigenvalues: +1 and -1. The hope is that incomplete factors will strike a good balance between computational cost and eigenvalue clustering.

It remains possible to use [MINRES](http://www.stanford.edu/group/SOL/software/minres.html) or [SYMMLQ](http://www.stanford.edu/group/SOL/software/symmlq.html) on the symmetrically-preconditioned system

|**D**|<sup>-1/2</sup> **L**<sup>-1</sup> **A**  **L**<sup>-T</sup> |**D**|<sup>-1/2</sup> **x** = |**D**|<sup>-1/2</sup> **L**<sup>-1</sup> **b**.

## Matlab Interface

The original Matlab interface has been fixed and modernized. Here is an example call:

    n = 6; m = 4; E = rand(m,n);
    A = [(n+m+1)*eye(n) E' ; E -(n+m+1)*eye(m)];
    Adiag = full(diag(A)); lA = sparse(tril(A,-1)); p=1;
    [L, Ldiag] = icfmex(lA, Adiag, p);
    L = L + speye(size(L,1));

## Notes

The original ICFS README is in `README.orig`.
