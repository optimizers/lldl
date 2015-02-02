# Julia interface to LLDL.

# liblldl.dylib should be in you LD_LIBRARY_PATH.
lldl_lib = "liblldl.dylib";
macro lldl_call(args...)
  quote
    ccall((:dicfs_, $lldl_lib), Void, $(args...))
  end
end

function lldl_base(K :: SparseMatrixCSC, p :: Int)
  n = size(K, 1);
  adiag = diag(K);
  A = sparse(tril(K, -1));
  nnzA = nnz(A);
  d = zeros(Float64, n);
  l = zeros(Float64, nnzA + n * p);
  lcolptr = zeros(Int32, n + 1);
  lrowind = zeros(Int32, nnzA + n * p);
  alpha = zeros(Float64, 1);  # May be modified by Fortran.

  # Temporary work arrays.
  iwa = zeros(Int32, 3 * n, 1);
  wa1 = zeros(Float64, n, 1);
  wa2 = zeros(Float64, n, 1);

  @lldl_call((Ptr{Int32},   Ptr{Int32},
              Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32},
              Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32},
              Ptr{Int32},   Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}),
              &n, &nnzA,
              A.nzval, adiag, convert(Array{Int32, 1}, A.colptr), convert(Array{Int32, 1}, A.rowval),
              l, d, lcolptr, lrowind,
              &p, alpha, iwa, wa1, wa2);

  return (l, lcolptr, lrowind, d, alpha[1]);
end

function lldl(K :: SparseMatrixCSC, p :: Int)
  n = size(K, 1);
  (l, lcolptr, lrowind, d, alpha) = lldl_base(K, p);
  L = SparseMatrixCSC(n, n, lcolptr, lrowind, l);
  return (L, d, alpha);
end

using LinearOperators

function lldl_op(K :: SparseMatrixCSC, p :: Int)
  n = size(K, 1);
  (Lmat, d, alpha) = lldl(K, p);
  D = opDiagonal(1./abs(d));
  L = opInverse(Lmat + speye(n));
  return (L' * D * L, alpha);
end
