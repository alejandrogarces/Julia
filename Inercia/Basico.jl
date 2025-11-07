using LinearAlgebra
using Plots
# Estado estacionario
pm = 8
x = 0.1
ξ = 1
H = 2
w_base = 2*pi*60
θ = asin(x*pm);
# Simulacion base
Δp = 0.1
np = 10000
dt = 0.0001
gr_ws = zeros(np)
gr_th = zeros(np)
gr_ws[1] = 1
gr_th[1] = θ
for k = 2:np
    gr_th[k] = gr_th[k-1] + dt*w_base*(gr_ws[k-1]-1)
    gr_ws[k] = gr_ws[k-1] + dt*(pm+Δp-1/x*sin(gr_th[k-1])-ξ*(gr_ws[k-1]-1))/(2*H)
end
t = (1:np)*dt
p1 = plot(t,gr_ws,label="ω")

# Eigenvalues
A = [-ξ/(2*H)  -cos(θ)/(2*H*x); w_base 0]
L = eigen(A)
p2 = plot(real.(L.values),imag.(L.values), seriestype=:scatter, label="λ")
p3 = plot(p1,p2,layout=(1,2),size = (614, 400),background_color = :transparent)
savefig("CasoBase.png")
display(p3)
# Cambio en la inercia
H = H + 3

# Transitorio
for k = 2:np
    gr_th[k] = gr_th[k-1] + dt*w_base*(gr_ws[k-1]-1)
    gr_ws[k] = gr_ws[k-1] + dt*(pm+Δp-1/x*sin(gr_th[k-1])-ξ*(gr_ws[k-1]-1))/(2*H)
end
p1 = plot!(p1,t,gr_ws,label="ω(+H)")

# Eigenvalues
A = [-ξ/(2*H)  -cos(θ)/(2*H*x); w_base 0]
L = eigen(A)
p2 = plot!(p2,real.(L.values),imag.(L.values), seriestype=:scatter, label="λ(+H)")
p3 = plot(p1,p2,layout=(1,2),size = (614, 400),background_color = :transparent)
savefig("CambioInercia.png")
display(p3)
# Cambio en el damping
ξ = ξ*20

# Transitorio
for k = 2:np
    gr_th[k] = gr_th[k-1] + dt*w_base*(gr_ws[k-1]-1)
    gr_ws[k] = gr_ws[k-1] + dt*(pm+Δp-1/x*sin(gr_th[k-1])-ξ*(gr_ws[k-1]-1))/(2*H)
end
p1 = plot!(p1,t,gr_ws,label="ω(+H+D)")

# Eigenvalues
A = [-ξ/(2*H)  -cos(θ)/(2*H*x); w_base 0]
L = eigen(A)
p2 = plot!(p2,real.(L.values),imag.(L.values), seriestype=:scatter, label="λ(+H+D)")
p3 = plot(p1,p2,layout=(1,2),size = (614, 400),background_color = :transparent)
savefig("CambioInerciayDamping.png")
display(p3)
(+H+D+x)")
p3 = plot(p1,p2,layout=(1,2),size = (614, 400),background_color = :transparent)
savefig("CambioInerciayDampingeImpedancia.png")
display(p3)
