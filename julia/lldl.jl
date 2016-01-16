# Julia interface to LLDL.

# liblldl.dylib should be in you LD_LIBRARY_PATH.
lldl_lib = "liblldl.dylib";
macro lldl_call(args...)
  quote
    ccall((:dicfs_, $lldl_lib), Void, $(args...))
  end
end

function lldl_base(K :: SparseMatrixCSC, p :: Int; droptol :: Float64=0.0, minpiv :: Float64=0.0)
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

  @printf("droptol=%8.1e, minpiv=%8.1e\n", droptol, minpiv)

  @lldl_call((Ptr{Int32},   Ptr{Int32},
              Ptr{Float64}, Ptr{Float64}, Ptr{Int32},   Ptr{Int32},
              Ptr{Float64}, Ptr{Float64}, Ptr{Int32},   Ptr{Int32},
              Ptr{Int32},   Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
              Ptr{Int32},   Ptr{Float64}, Ptr{Float64}),
              &n, &nnzA,
              A.nzval, adiag, convert(Array{Int32, 1}, A.colptr), convert(Array{Int32, 1}, A.rowval),
              l, d, lcolptr, lrowind,
              &p, alpha, &droptol, &minpiv,
              iwa, wa1, wa2);

  return (l, lcolptr, lrowind, d, alpha[1]);
end

function lldl(K :: SparseMatrixCSC, p :: Int; droptol :: Float64=0.0, minpiv :: Float64=0.0)
  n = size(K, 1);
  (l, lcolptr, lrowind, d, alpha) = lldl_base(K, p, droptol=droptol, minpiv=minpiv);
  L = SparseMatrixCSC(n, n, lcolptr, lrowind, l);
  return (L, d, alpha);
end

using LinearOperators

function lldl_op(K :: SparseMatrixCSC, p :: Int; droptol :: Float64=0.0, minpiv :: Float64=0.0)
  n = size(K, 1);
  (Lmat, d, alpha) = lldl(K, p, droptol=droptol, minpiv=minpiv);
  D = opDiagonal(1./abs(d));
  L = opInverse(Lmat + speye(n));
  return (L' * D * L, alpha);
end
