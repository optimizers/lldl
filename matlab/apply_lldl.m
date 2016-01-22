function x = apply_lldl(LLDL, p, b)
  %
  % x = apply_lldl(LLDL, p, b)
  %
  % Solve the system Ax=b using the incomplete LDL factorization of A
  % represented by the Spot operator `LLDL`. The argument `p` represents
  % a permutation vector such that, in matrix form, P'AP = LDL.
  %

  x = LLDL * b(p);
  x(p) = x;
end
