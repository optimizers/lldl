load(['../tprobs/1138bus.mat']);
Adiag = full(diag(A)); n = size(A,1); lA = tril(A,-1); p=5;
[L, Ldiag] = icfmex(n, lA, Adiag, 5);
