# LLDL

## Overview

LLDL is a modification of [Lin and Moré's limited-memory Cholesky factorization](http://dx.doi.org/10.1137/S1064827597327334) code-named ICFS for symmetric positive definite matrices. LLDL implements a similar limited-memory scheme for symmetric indefinite matrices that possess a LDL<sup>T</sup> factorization, i.e., with D diagonal. Symmetric quasi-definite matrices fall into this category.

LLDL is applicable to symmetric indefinite matrices that are not
quasi definite, or more generally do not admit a LDL<sup>T</sup>
factorization. In this case, it computes an incomplete LDL<sup>T</sup>
factorization of a nearby matrix.

The main idea is that if L and D are the exact factors of A, the preconditioned matrix (L |D| L<sup>T</sup>)<sup>-1</sup> A has only two eigenvalues: +1 and -1. The hope is that incomplete factors will strike a good balance between computational cost and eigenvalue clustering.

It remains possible to use [MINRES](http://www.stanford.edu/group/SOL/software/minres.html) or [SYMMLQ](http://www.stanford.edu/group/SOL/software/symmlq.html) on the symmetrically-preconditioned system

|D|<sup>-1/2</sup> L<sup>-1</sup> A  L<sup>-T</sup> |D|<sup>-1/2</sup> x = |D|<sup>-1/2</sup> L<sup>-1</sup> b.

## Installing

Just do it ;)

    ./install

You can use non-default compilers and compiler flags by passing options to the `install` script. For instance

    CC=clang FC=ifort CFLAGS='-g' FFLAGS='-FI' ./install

Here is the complete list:

    $ ./install --help
    install [--help] [--skip-matlab-check]
    The following environment variables may be set to influence
    the behavior of the install procedure:
      CC           C compiler executable
      FC           Fortran compiler executable
      CFLAGS       C compiler flags
      FFLAGS       Fortran compiler flags
      MATLABDIR    Path to Matlab's bin and extern subdirs

## Matlab Interface

You can specify the location of Matlab's `bin` and `extern` subdirectories using the `MATLABDIR` variable:

    MATLABDIR=/Applications/Matlab/MATLAB_R2012a.app ./install

The original Matlab interface has been fixed and modernized. The Mathworks [only support version 4.3 gfortran](http://www.mathworks.com/support/compilers/R2013a/index.html?sec=maci64) on OSX and Linux. On Linux, `gcc-4.3` should be found in your package manager. On OSX, I recommend using [Homebrew](http://mxcl.github.io/homebrew). Once Homebrew is installed, `gcc-4.3`, including the Fortran compiler, may be installed using

    brew tap homebrew/versions
    brew install gcc43 --enable-fortran --enable-profiled-build

The compiler executables installed by the above commands are those used by default in LLDL's `install` script. If you know that your compilers will produce valid MEX files (it is the case for `gcc-4.2` and `gfortran-4.2` on OSX 10.6.8), then:

    ./install --skip-matlab-check

Here is an example Matlab session:

    n = 6; m = 4; E = rand(m,n);
    A = [(n+m+1)*eye(n) E' ; E -(n+m+1)*eye(m)];
    Adiag = full(diag(A)); lA = sparse(tril(A,-1)); p=1;
    [L, D, shift] = lldl(lA, Adiag, p);
    L = L + speye(size(L));  % Diagonal was left out of L.

In practice it is more programmatically convenient to use the incomplete factors L and D as implicitly defining the preconditioner by bundling them into an abstract [SPOT](https://github.com/mpf/spot) operator:

    LLDL = opLLDL(A, p);
    x = LLDL * y;  % Solves the system L * |D| * L' x = y;

## Python Interface

A Python interface is included as part of [NLPy](https://github.com/dpo/nlpy).

## Trouble / Questions / Bugs

Search existing [issues](https://github.com/optimizers/lldl/issues). If that does not answer your question, please open a [new issue](https://github.com/optimizers/lldl/issues/new).

## To Do List

This is a list of improvements that should be added to LLDL to improve its efficiency, particularly with Matlab.

- [ ] Implicit permutation: instead of having to explicity apply a reordering, the user supplies a permutation vector. The incomplete factorization proceeds in the order given by this permutation vector.
- [x] In a typical Matlab application, only the lower triangle of the matrix should be kept in memory and should be wrapped into a SPOT operator to compute matrix-vector products and take the permutation into account implicitly.
- [x] Reinstate the modification of the diagonal to obtain a SQD H-matrix. A first idea is A(i,i) += α if A(i,i) > 0 and A(i,i) -= α if A(i,i) < 0, where α > 0.

## References

* [C.-J. Lin and J. J. Moré. Incomplete Cholesky factorizations with limited memory. *SIAM Journal on Scientific Computing*, 21(1):24--45, 1999](http://dx.doi.org/10.1137/S1064827597327334).

## Notes

The original ICFS README is in `README.orig`.
