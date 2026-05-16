# Modelo de transportes
# alejandro.garces@utp.edu.co

using LinearAlgebra
using JuMP
using HiGHS

# Parametros
gmax = [150;  0;360;  0;  0;600]
d = [ 80;240; 40;160;240;  0]
lineas = [(1,2);(1,4);(1,5);(2,3);(2,4);(2,6);(3,5);(4,6)]
c = [40;60;20;20;40;30;20;30]
x = [0.4;0.6;0.2;0.2;0.4;0.3;0.2;0.3]
fmax = [100;80;100;100;100;100;100;100]
z0 = [1;1;1;1;1;0;1;0]

# triplicar las lineas para adminir mas líneas por corredor
lineas=repeat(lineas,3)
c = repeat(c,3)
x = repeat(x,3)
fmax = repeat(fmax,3)
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

# Modelo de transportes
ModeloTransportes = Model(HiGHS.Optimizer)
@variable(ModeloTransportes, f[1:nlin])
@variable(ModeloTransportes, g[1:nnod]>=0)
@variable(ModeloTransportes, z[1:nlin]>=0, Bin)
@objective(ModeloTransportes, Min, c'*z)
@constraint(ModeloTransportes,S*f+g==d)
@constraint(ModeloTransportes, f .<= z.*fmax)
@constraint(ModeloTransportes, f .>= -z.*fmax)
@constraint(ModeloTransportes, g .<= gmax)
@constraint(ModeloTransportes, z[1:ncor] == z0)

optimize!(ModeloTransportes)

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
