# Planeacion de la generacion
#  alejandro.garces@utp.edu.co
using LinearAlgebra
using JuMP
using HiGHS
using DataFrames
# Parametros
d = 9
data = DataFrame()
data[!,:Marca] = ["LG","GE","MG","HW"]
data[!,:α] = [100,160,200,300]
data[!,:β] = [51, 45, 43, 40]
data[!,:pmax] = [5,5,5,5]
println(data)
n = length(data[:,:Marca])
# Modelo de optimizacion
PlanGen = Model(HiGHS.Optimizer)
@variable(PlanGen, p[1:n]>=0)
@variable(PlanGen, z[1:n],Bin)
@objective(PlanGen,Min, sum(data[:,:α].*z[:])+sum(data[:,:β].*p[:]))
@constraint(PlanGen,sum(p)==d)
@constraint(PlanGen,p.<=data[:,:pmax].*z[:])
optimize!(PlanGen)

# Mostrar resultados
data[:,:pgen] = value.(p)
data[:,:Comprar] = value.(z)
println(data)
println(objective_value(PlanGen))
