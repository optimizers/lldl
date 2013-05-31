# LLDL

## Overview

LLDL is a modification of [Lin and Moré's limited-memory Cholesky factorization](http://dx.doi.org/10.1137/S1064827597327334) code-named ICFS for symmetric positive definite matrices. LLDL implements a similar limited-memory scheme for symmetric indefinite matrices that possess a LDL<sup>T</sup> factorization, i.e., with D diagonal. Symmetric quasi-definite matrices fall into this category.

The main idea is that if L and D are the exact factors of A, the preconditioned matrix (L |D| L<sup>T</sup>)<sup>-1</sup> A has only two eigenvalues: +1 and -1. The hope is that incomplete factors will strike a good balance between computational cost and eigenvalue clustering.

It remains possible to use [MINRES](http://www.stanford.edu/group/SOL/software/minres.html) or [SYMMLQ](http://www.stanford.edu/group/SOL/software/symmlq.html) on the symmetrically-preconditioned system

|D|<sup>-1/2</sup> L<sup>-1</sup> A  L<sup>-T</sup> |D|<sup>-1/2</sup> x = |D|<sup>-1/2</sup> L<sup>-1</sup> b.

## Installing

Just do it ;)

    ./install

You can use non-default compilers and compiler flags by passing options to the `install` script. For instance

    CC=clang FC=ifort CFLAGS='-g' FFLAGS='-FI' ./install

## Matlab Interface

The original Matlab interface has been fixed and modernized. The Mathworks [only support version 4.3 gfortran](http://www.mathworks.com/support/compilers/R2013a/index.html?sec=maci64) on OSX and Linux. On Linux, `gcc-4.3` should be found in your package manager. On OSX, I recommend using [Homebrew](http://mxcl.github.io/homebrew). Once Homebrew is installed, `gcc-4.3`, including the Fortran compiler, may be installed using

    brew tap homebrew/versions
    brew install gcc43 --enable-fortran --enable-profiled-build

The compiler executables installed by the above commands are those used by default in LLDL's `install` script.

Here is an example call:

    n = 6; m = 4; E = rand(m,n);
    A = [(n+m+1)*eye(n) E' ; E -(n+m+1)*eye(m)];
    Adiag = full(diag(A)); lA = sparse(tril(A,-1)); p=1;
    [L, Ldiag] = icfmex(lA, Adiag, p);
    L = L + speye(size(L,1));

## Python Interface

A Python interface is included as part of [NLPy](https://github.com/dpo/nlpy).

## Notes

The original ICFS README is in `README.orig`.
