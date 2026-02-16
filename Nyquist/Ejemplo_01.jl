using LinearAlgebra
using Plots
using Pkg
using Polynomials
theme(:dark)
np = 200
w = (0:np)*0.01
s = w*1im
p1 = 0.1+1im
p2 = conj(p1)
p3 = 0.1
K = 1
G(s) = K/((s+p1)*(s+p2)*(s+p3))
zp = G.(s)
zn = G.(-s)
z = [zp;zn[np:-1:1]]
pol1 = [p1*p2*p3, p1*p2+p1*p3+p2*p3,p1+p2+p3,1]
pol2 = [p1*p2*p3+K, p1*p2+p1*p3+p2*p3,p1+p2+p3,1]
println("Denominador =", Polynomial(pol1))
println("Lazo cerrado = ", maximum(real.(roots(Polynomial(pol2)))))
n = length(z)
scatter([-1],[0], label="", xlimits=[-5.5,10.5],ylimits=[-5.5,5.5])
@gif for k in 1:n
    plot!(real.(z[1:k]),imag.(z[1:k]),label="",color=:orange, linewidth=2)
end
