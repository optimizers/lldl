beep off;
load(['../tprobs/1138bus.mat']);
Adiag = full(diag(A)); lA = tril(A,-1); p=5;
[L, Ldiag] = icfmex(lA, Adiag, 5);
