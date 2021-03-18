# Julia interface to LLDL.
using LinearOperators

# liblldl.dylib should be in you LD_LIBRARY_PATH.
const liblldl = "liblldl.dylib"
macro lldl_call(args...)
  args = map(esc, args)
  quote
    ccall((:dicfs_, liblldl), Nothing, $(args...))
  end
end

function lldl_base(K::SparseMatrixCSC, p::Int; droptol::Float64=0.0, minpiv::Float64=0.0)
  n = size(K, 1)
  adiag = convert(Vector{Float64}, diag(K))
  A = sparse(tril(K, -1))
  nnzA = nnz(A)
  d = zeros(n)
  l = zeros(nnzA + n * p)
  lcolptr = zeros(Int32, n + 1)
  lrowind = zeros(Int32, nnzA + n * p)
  alpha = zeros(1)  # May be modified by Fortran.

  # Temporary work arrays.
  iwa = zeros(Int32, 3 * n)
  wa1 = zeros(n)
  wa2 = zeros(n)

  @lldl_call((Ref{Int32},   Ref{Int32},
              Ref{Float64}, Ref{Float64}, Ref{Int32},   Ref{Int32},
              Ref{Float64}, Ref{Float64}, Ref{Int32},   Ref{Int32},
              Ref{Int32},   Ref{Float64}, Ref{Float64}, Ref{Float64},
              Ref{Int32},   Ref{Float64}, Ref{Float64}),
              n, nnzA,
              A.nzval, adiag, convert(Vector{Int32}, A.colptr), convert(Vector{Int32}, A.rowval),
              l, d, lcolptr, lrowind,
              p, alpha, droptol, minpiv,
              iwa, wa1, wa2)

  return (l, lcolptr, lrowind, d, alpha[1])
end

function lldl(K::SparseMatrixCSC, p::Int; droptol::Float64=0.0, minpiv::Float64=0.0)
  n = size(K, 1)
  (l, lcolptr, lrowind, d, alpha) = lldl_base(K, p, droptol=droptol, minpiv=minpiv)
  L = SparseMatrixCSC(n, n, lcolptr, lrowind, l)
  return (L, d, alpha)
end

function lldl_op(K::SparseMatrixCSC, p::Int; droptol::Float64=0.0, minpiv::Float64=0.0)
  n = size(K, 1)
  (Lmat, d, alpha) = lldl(K, p, droptol=droptol, minpiv=minpiv)
  D = opDiagonal(1 ./ abs.(d))
  L = opInverse(Lmat + spdiagm(0 => ones(n)))
  return (L' * D * L, alpha)
end
