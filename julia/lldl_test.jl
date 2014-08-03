include("lldl.jl")

n = 10;
m = 6;
E = sprand(m, n, .2);
K = [(n+m+1)*speye(n) E' ; E  -(n+m+1)*speye(m)];

@printf("p  ‖LDL'-K‖  κ(P*K)\n");
for p = 1 : 5
  (L, d, alpha) = lldl(K, p); L = L + speye(n+m);
  (P, alpha) = lldl_op(K, p);
  @printf("%1d  %8.2e  %8.2e\n", p, norm(L*diagm(d)*L' - K), cond(full(P * K)));
end

