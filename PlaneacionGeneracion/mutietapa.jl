# Planeacion de la generacion multietapa
#  alejandro.garces@utp.edu.co
using LinearAlgebra
using JuMP
using HiGHS
using DataFrames
# Parametros
d = [3,5,7,8,8,9]
data = DataFrame()
data[!,:Marca] = ["LG","GE","MG","HW"]
data[!,:α] = [100,160,200,300]
data[!,:β] = [51, 45, 43, 40]
data[!,:pmax] = [5,5,5,5]
println(data)
n = length(data[:,:Marca])
nt = length(d)
# Modelo de optimizacion
PlanGen = Model(HiGHS.Optimizer)
@variable(PlanGen, p[1:n,1:nt]>=0)
@variable(PlanGen, z[1:n,1:nt],Bin) # indica si esta on 
@variable(PlanGen, y[1:n,1:nt],Bin) # indica cuando se hizo la inversion
@objective(PlanGen,Min, sum(y'*data[:,:α])+sum(p'*data[:,:β]))
for t = 1:nt
    @constraint(PlanGen,sum(p[:,t])==d[t])
    @constraint(PlanGen,p[:,t].<=data[:,:pmax].*z[:,t])
end
@constraint(PlanGen,y[:,1] == z[:,1])
for t = 2:nt
    @constraint(PlanGen,y[:,t] == z[:,t]-z[:,t-1])
end
optimize!(PlanGen)

# Mostrar resultados
for t = 1:nt
  data[:,"pgenT"*string(t)] = value.(p[:,t])
  #data[:,:Comprar] = value.(z)
end
println(data)
println(objective_value(PlanGen))
