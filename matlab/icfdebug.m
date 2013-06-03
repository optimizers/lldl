beep off;
format short e
randseed;

%load(['../tprobs/1138bus.mat']);
n = 6; m = 4; E = rand(m,n);
A = [(n+m+1)*eye(n) E' ; E -(n+m+1)*eye(m)];
%A = 2*speye(n+m);
%A = [1 0 1 ; 0 1 1 ; 1 1 -1 ];
Adiag = full(diag(A)); lA = sparse(tril(A,-1)); p=0;  % Compute exact factor.

% [L, D] = icfmex(lA, Adiag, p);
% L = L + speye(size(L,1));       % L was the strict lower triangle only.

% A - full(L) * diag(D) * full(L')

for p = 0 : 4
  [L, D] = icfmex(lA, Adiag, p);
  L = L + speye(size(L,1));       % L was the strict lower triangle only.
  full(L)
end
