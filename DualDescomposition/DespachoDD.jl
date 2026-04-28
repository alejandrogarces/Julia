"""
    Metodo de descomposicion dual para
    Por: Alejandro Garces Ruiz
    28/04/2026
"""

using LinearAlgebra
using DataFrames

# Parametros
df = DataFrame()
df.a = [0.30494, 0.21174, 0.07092, 0.05606, 0.03598, 0.04222]
df.b = [38.5390, 46.1591, 38.3055, 40.3965, 38.2704, 36.3278]
df.pmax = [100, 100, 200, 200, 300, 500]
demanda = 1200
α = 0.025
println(df)

# Nodos esclavo
function DespachoUnidadTermica(λ,k)
    p = -(df.b[k]+λ)/df.a[k]
    p = minimum([p,df.pmax[k]])
    p = maximum([p,0])
    return p
end


# Nodo Maestro
n = length(df.pmax)
pdesp = zeros(n)
global λ = 0
for i = 1:10
    for k = 1:n
        pdesp[k] = DespachoUnidadTermica(λ,k)        
    end
    println("Iteracion ",i, " p = ", round.(pdesp,digits=2))
    global λ = λ + α*(sum(pdesp)-demanda)
end
