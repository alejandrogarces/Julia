# Generacion de escenarios usando Kmeans
#   alejandro.garces@utp.edu.co

using CSV 
using DataFrames 
using Plots
using LinearAlgebra
using Clustering
theme(:dark)

Solar = DataFrame(CSV.File("Texas_Radiation.csv"))
ncol = 365
nfil = div(length(Solar[:,"Power(MW)"]),ncol)
Data = reshape(Solar[:,"Power(MW)"], nfil, ncol)

#plot(Solar[:,"Power(MW)"])
plot(Data,label="",linecolor=:gray)

# Uso de kmeans
Escenarios = kmeans(Data,3)
Prob = round.(Escenarios.counts/ncol*100, digits=2)
println("Probabilidades por escenario", Prob)
plot!(Escenarios.centers, label= string.(Prob'))
