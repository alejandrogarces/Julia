# Ejemplo simple de ajuste de curvas usando optimizacion\
#   alejandro.garces@utp.edu.co
using Plots
using LinearAlgebra
using JuMP
using Ipopt
# Datos
n = 20
x = [k for k = 1:n];
y = 2*x .+ 0.2*randn(n) .+ n*ones(1)
y[Integer.([0.5,0.7,0.8].*n)] += 4*n*rand(3)  # datos espurio
y = round.(y,digits=3)
# Minimos cuadrados convencional
modelo1 = Model(Ipopt.Optimizer)
set_silent(modelo1)
@variable(modelo1, a1)
@variable(modelo1, b1)
@objective(modelo1,Min,sum(1/2*(y.-a1*x .-b1).^2))
optimize!(modelo1)
a1 = value(a1)
b1 = value(b1)
# norma-1
modelo2 = Model(Ipopt.Optimizer)
set_silent(modelo2)
@variable(modelo2, a2)
@variable(modelo2, b2)
@objective(modelo2,Min,sum(abs.(y.-a2*x.-b2)))
optimize!(modelo2)
a2 = value(a2)
b2 = value(b2)
#  Pseudo-Huber loss function 
delta = 0.01
modelo3 = Model(Ipopt.Optimizer)
set_silent(modelo3)
@variable(modelo3, a3)
@variable(modelo3, b3)
@objective(modelo3,Min,sum([delta^2*sqrt((1+(y[k]-a3*x[k]-b3)^2/delta^2)-1) for k = 1:n]))
optimize!(modelo3)
a3 = value(a3)
b3 = value(b3)
# Graficos
plt = scatter(x,y, label="datos")
plt = plot!(x,a1*x.+b1,label="Norma-2 a="*string(round(a1,digits=3)))
plt = plot!(x,a2*x.+b2,label="Norma-1 a="*string(round(a2,digits=3)))
plt = plot!(x,a3*x.+b3,label="Huber-M a="*string(round(a3,digits=3)))
