"""
**Ejemplo simple de planeación de la distribución**

    alejandro.garces@utp.edu.co
"""

# Instalar las librerias
## Instalar paquekes
paquetes = ["Plots", "LinearAlgebra", "JuMP", "DataFrames", "HiGHS"]
using Pkg
for paquete in paquetes
    if !haskey(Pkg.project().dependencies, paquete)
        @info "Instalando $paquete..."
        Pkg.add(paquete)
    end
end
for paquete in paquetes
	 @eval using $(Symbol(paquete))
end
theme(:dark)

# Datos de los nodos
Nodes = DataFrame()
Nodes.x = [8, 15, 17, 10, 7, 8, 17, 16, 10, 20, 21, 12, 2, 1, 10, 11, 13, 4, 13, 4, 10, 10, 15, 4]
Nodes.y = [-4, -3, -6, -10, -11, -16, -15, -18, -7, -10, -13, -5, -13, -7, -13, -15, -7, -4, -14, -16, -2, -17, -10, -10]
Nodes.On   = zeros(24)
Nodes.Dem = [2.54,0.50,1.66,0.29,0.18,0.75,2.87,0.46,0.80,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
Nodes.Sub = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1]
Nodes.Cost = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,5,5,5]
println(Nodes)

# Datos de los tramos de línea
Lines = DataFrame()
Lines.Bus1 = [0,0,0,0,0,20,1,1,8,2,2,2,3,3,3,4,4,5,5,5,6,6,6,6,6,7,21,17,13,23,12,12,14,14,16,18,9,22] .+ 1
Lines.Bus2 = [20,8,11,17,4,1,2,11,3,16,9,22,14,16,6,5,23,21,19,15,7,9,22,10,18,21,15,13,23,12,19,14,15,18,9,10,22,10] .+ 1
Lines.On =   zeros(38)
Lines.Smax = 4*ones(38)
Lines.Cost = 0.1*ones(38)
println(Lines)

# Grafico
function GraficarRed(Nodes,Lines)
  NumN = length(Nodes.x)
  NumL = length(Lines.Bus1)
  NumS = sum(Nodes.Sub)
  plt = scatter(Nodes.x[NumN-NumS+1:NumN],Nodes.y[NumN-NumS+1:NumN],legend=false,markersize = 16,color=:black, size=(600,600), axis=false)
  plt = scatter!(Nodes.x,Nodes.y,legend=false,markersize = 12,color=:gray)
  for k = 1:NumL
      n1 = Lines.Bus1[k]
      n2 = Lines.Bus2[k]
      if Lines.On[k] != 0
        plt = plot!(Nodes.x[[n1,n2]],Nodes.y[[n1,n2]], linewidth=2, linecolor = :orange)
      else
        plt = plot!(Nodes.x[[n1,n2]],Nodes.y[[n1,n2]], linestyle = :dash, linecolor=:gray)
      end
  end
  for k = 1:NumN
      if Nodes.Dem[k] > 0 # nodos de carga
        plt = scatter!([Nodes.x[k]],[Nodes.y[k]], legend=false, markersize = 16, color=:gray)
      end
      if Nodes.On[k] != 0
        plt = scatter!([Nodes.x[k]],[Nodes.y[k]], legend=false, markersize = 12, color=:orange)
      end
      plt = annotate!(Nodes.x[k],Nodes.y[k],(k,:black,8,:center))
  end
return plt
end
GraficarRed(Nodes,Lines)

# Matriz de incidencia
function MatrizIncidencia(Lines,Nodes)
    NumN = length(Nodes.x)
    NumL = length(Lines.Bus1)
    A = zeros(NumN,NumL)
    for k = 1:NumL
        n1 = Lines.Bus1[k]
        n2 = Lines.Bus2[k]
        A[n1,k] = 1
        A[n2,k] = -1
    end
return A
end
A=MatrizIncidencia(Lines,Nodes)

# Modelo de transportes
NumS = sum(Nodes.Sub) # Numero de subestaciones (las NumS últimas)
NumD = sum(Nodes.Dem.>0) # Numero de nodos con demanda
NumN = length(Nodes.x)
NumL = length(Lines.Bus1)
Pt = sum(Nodes.Dem)
PlanDist = Model(HiGHS.Optimizer)
@variable(PlanDist,g[1:NumN]>=0)
@variable(PlanDist,p[1:NumL])
@variable(PlanDist,zN[1:NumN],Bin)
@variable(PlanDist,zL[1:NumL],Bin)
@variable(PlanDist,f[1:NumL]) # flujo ficticio
@objective(PlanDist,Min,sum(Nodes.Cost.*zN)+sum(Lines.Cost.*zL))
@constraint(PlanDist,g .<= Nodes.Sub.*zN*Pt)
@constraint(PlanDist,p .<=  Lines.Smax.*zL)
@constraint(PlanDist,p .>= -Lines.Smax.*zL)
@constraint(PlanDist,g-Nodes.Dem == A*p)
@constraint(PlanDist,sum(zL) == sum(zN)-1)
zS = [zeros(NumN-NumS,1);zN[(NumN-NumS+1):NumN]] # subestaciones
@constraint(PlanDist,sum(zS) == 1) # Solo una sola subestacion
#@constraint(PlanDist,NumN*zS[:]-zN[:] = A*f)
@constraint(PlanDist,NumN*zS .- zN .>= A*f)
@constraint(PlanDist,f .>= -(NumN-1).*zL)
@constraint(PlanDist,f .<=  (NumN-1).*zL)
optimize!(PlanDist)

ResN = Nodes
ResL = Lines
ResN.On = round.(value.(zN))
ResL.On = round.(value.(zL))
GraficarRed(ResN,ResL)
