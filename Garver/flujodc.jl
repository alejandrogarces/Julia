# Modelo de transportes
# alejandro.garces@utp.edu.co

using LinearAlgebra
using JuMP
using HiGHS

# Parametros
gmax = [150;  0;360;  0;  0;600]
#gmax = [50; 0; 165; 0; 0; 600]
d = [ 80;240; 40;160;240;  0]
lineas = [(1,2);(1,4);(1,5);(2,3);(2,4);(2,6);(3,5);(4,6)]
c = [40;60;20;20;40;30;20;30]
x = [0.4;0.6;0.2;0.2;0.4;0.3;0.2;0.3]
fmax = [100;80;100;100;100;100;100;100]
z0 = [1;1;1;1;1;0;1;0]

# triplicar las lineas para adminir mas líneas por corredor
nl = 10
lineas=repeat(lineas,nl)
c = repeat(c,nl)
x = repeat(x,nl)
fmax = repeat(fmax,nl)
ncor = length(z0)    
nlin = length(c)
nnod = length(d)
c[1:ncor] .= 0. # los corredores existentes tienen costo cero

# Matriz de incidencia
S = zeros(nnod,nlin)
for (lin,conexion) in enumerate(lineas)
    k = Int(conexion[1])
    m = Int(conexion[2])
    S[k,lin] = 1
    S[m,lin] = -1
end

# Modelo de con flujo DC
η = pi
ModeloDC = Model(HiGHS.Optimizer)
@variable(ModeloDC, f[1:nlin])
@variable(ModeloDC, g[1:nnod]>=0)
@variable(ModeloDC, θ[1:nnod]>=0)
@variable(ModeloDC, z[1:nlin], Bin)
@objective(ModeloDC, Min, c'*z)
@constraint(ModeloDC,S*f+g==d)
@constraint(ModeloDC, f .<= z.*fmax)
@constraint(ModeloDC, f .>= -z.*fmax)
@constraint(ModeloDC, g .<= gmax)
@constraint(ModeloDC, z[1:ncor] == z0)
@constraint(ModeloDC,-z*η <= x.*f/100)
@constraint(ModeloDC, x.*f/100 .<= z*η )
for (k,conexion) in enumerate(lineas)
    i = Int(conexion[1])
    j = Int(conexion[2])
    @constraint(ModeloDC, -(1-z[k])*η <= x[k].*f[k]/100 - (θ[i]-θ[j]) )
    @constraint(ModeloDC, x[k].*f[k]/100 - (θ[i]-θ[j]) <= (1-z[k])*η )
end
optimize!(ModeloDC)

println("Costo:", value(c'*z))
println("Generación")
for i = 1:nnod
    if value(g[i])>=0.1
        println("g",i,"=",value(g[i]))
    end 
end
println("Nuevas líneas")
for k = 1:nlin
    if (value(z[k])>=0.1)&(c[k]>0)
        println(lineas[k],"=",round(value(z[k])))
    end 
end
