include("lldl.jl")

n = 6;
m = 4;
E = sprand(m, n, .2);
K = [(n+m+1)*speye(n) E' ; E  -(n+m+1)*speye(m)];
#=K = [(n+m+1)*eye(n) E' ; E  -(n+m+1)*zeros(m,m)];=#
p = 1;

(L, d, alpha) = lldl_mat(K, p);
