include("lldl.jl")
using Printf

n = 10
m = 6
E = sprand(m, n, .2)
In = (n + m + 1) * spdiagm(0 => ones(n))
Im = (n + m + 1) * spdiagm(0 => ones(m))
K = [In E' ; E  -Im]

@printf("p  ‖LDL'-K‖  κ(P*K)\n")
for p = 1 : 5
  (L, d, alpha) = lldl(K, p)
  L = L + I
  (P, alpha) = lldl_op(K, p)
  @printf("%1d  %8.2e  %8.2e\n", p, norm(L*diagm(d)*L' - K), cond(Matrix(P * K)))
end

