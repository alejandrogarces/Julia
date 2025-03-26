"""
  Flujo de carga óptimo para un sistema de dos nodos
  en donde la potencia en el nodo de recibo es -3.5-1.5j
  y la impedancia de las líneas es de 1/(0.010+0.023)
  Comparacion entre diferentes formas de convexificar la funcion
"""
using LinearAlgebra
using JuMP
using Ipopt
using DataFrames
r = 0.010
x = 0.023
yL = 1/(r+x*1im)
g = real(yL)
b = imag(yL)
p_esp_m = -3.5
q_esp_m = -1.5
results = DataFrame()
results[!,"Modelo"] = ["Completo";"Wirtinger";"Secuencial";"SOC";"SDP"]
results[!,"V1"] = [0.0;0;0;0;0]
results[!,"V2"] = [0.0;0;0;0;0]
results[!,"P1"] = [0.0;0;0;0;0]
results[!,"P2"] = [0.0;0;0;0;0]
results[!,"Q1"] = [0.0;0;0;0;0]
results[!,"Q2"] = [0.0;0;0;0;0]
results[!,"Perdidas"] = [0.0;0;0;0;0]

# Modelo no-lineal exacto
println("Modelo exacto")
model = Model(Ipopt.Optimizer)
@variable(model,vk_real)
@variable(model,vk_imag)
@variable(model,vm_real)
@variable(model,vm_imag)
@variable(model,pk)
@variable(model,qk)
@variable(model,pm)
@variable(model,qm)
vk = vk_real+vk_imag*1im
vm = vm_real+vm_imag*1im
vkm = vk-vm
@objective(model, Min, real(yL*vkm*conj(vkm)))
@NLconstraint(model, vk_real^2+vk_imag^2<=1.05^2)
@NLconstraint(model, vk_real^2+vk_imag^2>=0.95^2)
@constraint(model,vk_real>=0.95)
@constraint(model,vk_real<=1.05)
@NLconstraint(model, vm_real^2+vm_imag^2<=1.05^2)
@NLconstraint(model, vm_real^2+vm_imag^2>=0.95^2)
@constraint(model,vm_real>=0.95)
@constraint(model, conj(pk+qk*1im) ==  yL*conj(vk)*vkm)
@constraint(model, conj(pm+qm*1im) == -yL*conj(vm)*vkm)
@constraint(model, pm == p_esp_m)
@constraint(model, qm == q_esp_m)
set_silent(model)
optimize!(model)
println("\t",termination_status(model))
results."V1"[1] = sqrt(value(vk_real^2+vk_imag^2))
results."V2"[1] = sqrt(value(vm_real^2+vm_imag^2))
results."P1"[1] = value(pk)
results."P2"[1] = value(pm)
results."Q1"[1] = value(qk)
results."Q2"[1] = value(qm)
results."Perdidas"[1] = value(sum(pk+pm))

# Linealizacion
println("Linealizacion")
uk = 1+0im
um = 1+0im
sk_conj = yL*conj(uk)*(uk-um)
sm_conj = yL*conj(um)*(um-uk)
model2 = Model(Ipopt.Optimizer)
@variable(model2,vk_real)
@variable(model2,vk_imag)
@variable(model2,vm_real)
@variable(model2,vm_imag)
@variable(model2,pk)
@variable(model2,qk)
@variable(model2,pm)
@variable(model2,qm)
vk = vk_real+vk_imag*1im
vm = vm_real+vm_imag*1im
vkm = vk-vm
Δvk = vk-uk
Δvm = vm-um

@objective(model2, Min, real(yL*vkm*conj(vkm)))
@constraint(model2,vk_real>= 0.95)
@constraint(model2,vk_real<= 1.05)
@constraint(model2,vk_imag>=-0.05)
@constraint(model2,vk_imag<= 0.05)
@constraint(model2,vm_real>= 0.95)
@constraint(model2,vm_real<= 1.05)
@constraint(model2,vm_imag>=-0.05)
@constraint(model2,vm_imag<= 0.05)
@NLconstraint(model2, vk_real^2+vk_imag^2<=1.05^2)
@constraint(model2, qm == q_esp_m)
@constraint(model2, conj(pk+qk*1im)-sk_conj ==  yL*conj(uk)*(Δvk-Δvm) + yL*(uk-um)*conj(Δvk))
@constraint(model2, conj(pm+qm*1im)-sm_conj ==  yL*conj(um)*(Δvm-Δvk) + yL*(um-uk)*conj(Δvm))
@constraint(model2, pm == p_esp_m)
set_silent(model2)
optimize!(model2)
println("\t",termination_status(model2))
results."V1"[2] = sqrt(value(vk_real^2+vk_imag^2))
results."V2"[2] = sqrt(value(vm_real^2+vm_imag^2))
results."P1"[2] = value(pk)
results."P2"[2] = value(pm)
results."Q1"[2] = value(qk)
results."Q2"[2] = value(qm)
results."Perdidas"[2] = value(sum(pk+pm))

# Linealizacion segunda linealizacion
uk = value(vk)
um = value(vm)
sk_conj = yL*conj(uk)*(uk-um)
sm_conj = yL*conj(um)*(um-uk)
model2 = Model(Ipopt.Optimizer)
@variable(model2,vk_real)
@variable(model2,vk_imag)
@variable(model2,vm_real)
@variable(model2,vm_imag)
@variable(model2,pk)
@variable(model2,qk)
@variable(model2,pm)
@variable(model2,qm)
vk = vk_real+vk_imag*1im
vm = vm_real+vm_imag*1im
vkm = vk-vm
Δvk = vk-uk
Δvm = vm-um

@objective(model2, Min, real(yL*vkm*conj(vkm)))
@constraint(model2,vk_real>= 0.95)
@constraint(model2,vk_real<= 1.05)
@constraint(model2,vk_imag>=-0.05)
@constraint(model2,vk_imag<= 0.05)
@constraint(model2,vm_real>= 0.95)
@constraint(model2,vm_real<= 1.05)
@constraint(model2,vm_imag>=-0.05)
@constraint(model2,vm_imag<= 0.05)
@NLconstraint(model2, vk_real^2+vk_imag^2<=1.05^2)
@constraint(model2, qm == q_esp_m)
@constraint(model2, conj(pk+qk*1im)-sk_conj ==  yL*conj(uk)*(Δvk-Δvm) + yL*(uk-um)*conj(Δvk))
@constraint(model2, conj(pm+qm*1im)-sm_conj ==  yL*conj(um)*(Δvm-Δvk) + yL*(um-uk)*conj(Δvm))
@constraint(model2, pm == p_esp_m)
set_silent(model2)
optimize!(model2)
println("Modelo aproximado")
println("\t",termination_status(model2))
results."V1"[3] = sqrt(value(vk_real^2+vk_imag^2))
results."V2"[3] = sqrt(value(vm_real^2+vm_imag^2))
results."P1"[3] = value(pk)
results."P2"[3] = value(pm)
results."Q1"[3] = value(qk)
results."Q2"[3] = value(qm)
results."Perdidas"[3] = value(sum(pk+pm))

# Modelo SOC
modelsoc = Model(Ipopt.Optimizer)
@variable(modelsoc,uk)
@variable(modelsoc,um)
@variable(modelsoc,wkm_real)
@variable(modelsoc,wkm_imag)
@variable(modelsoc,wmk_real)
@variable(modelsoc,wmk_imag)
@variable(modelsoc,pk)
@variable(modelsoc,qk)
@variable(modelsoc,pm)
@variable(modelsoc,qm)
wkm = wkm_real + wkm_imag*1im
wmk = wmk_real + wmk_imag*1im

@objective(modelsoc, Min, real(yL*(uk+um-wkm-wmk)))
@constraint(modelsoc, uk >=0.9025)
@constraint(modelsoc, um >=0.9025)
@constraint(modelsoc, uk <=1.1025)
@constraint(modelsoc, um <=1.1025)
@constraint(modelsoc, um <=1.1025)
@constraint(modelsoc, wkm == conj(wmk))
@NLconstraint(modelsoc, (4*(wkm_real^2+wkm_imag^2)+(uk-um)^2) <= (uk+um)^2)
@constraint(modelsoc, conj(pk+qk*1im) == yL*(uk-wkm))
@constraint(modelsoc, conj(pm+qm*1im) == yL*(um-wmk))
@constraint(modelsoc, pm == p_esp_m)
@constraint(modelsoc, qm == q_esp_m)
set_silent(modelsoc)
optimize!(modelsoc)
println("Modelo SOC")
println("\t",termination_status(modelsoc))
results."V1"[4] = sqrt(value(uk))
results."V2"[4] = sqrt(value(um))
results."P1"[4] = value(pk)
results."P2"[4] = value(pm)
results."Q1"[4] = value(qk)
results."Q2"[4] = value(qm)
results."Perdidas"[4] = value(sum(pk+pm))

# Modelo SDP
import SCS
modelo_sdp = Model(SCS.Optimizer)
Y = [yL -yL;-yL yL+0.001+0.001im]  # evitar que la matriz sea singular
G = real(Y)
B = imag(Y)
iG = pinv(G)
iB = pinv(B)
@variable(modelo_sdp, X[1:2, 1:2], PSD)
@variable(modelo_sdp, M[1:2, 1:2])
@variable(modelo_sdp, N[1:2, 1:2])
@variable(modelo_sdp, P[1:2])
@variable(modelo_sdp, Q[1:2])
@constraint(modelo_sdp,  M  + N == (iB*G+iG*B)*X)
@constraint(modelo_sdp,  M' - N == X*(iB*G)' - (iG*B)*X)
@constraint(modelo_sdp,   P == diag(B*M))
@constraint(modelo_sdp,  -Q == diag(G*N))
@constraint(modelo_sdp, P[2] == p_esp_m)
@constraint(modelo_sdp, Q[2] == q_esp_m)
@constraint(modelo_sdp, X.>= 0.9025)
@constraint(modelo_sdp, X.<= 1.1025)
@objective(modelo_sdp, Min, tr(G*X))
println("Modelo SDP")
set_silent(modelo_sdp)
optimize!(modelo_sdp)
println("\t",termination_status(modelo_sdp))
L = eigvals(value(X))
M = eigvecs(value(X))
V = abs.(sqrt(L[2])*M[:,2])
results."V1"[5] = V[1]
results."V2"[5] = V[2]
results."P1"[5] = value(P[1])
results."P2"[5] = value(P[2])
results."Q1"[5] = value(Q[1])
results."Q2"[5] = value(Q[2])
results."Perdidas"[5] = value(sum(P))
println(results)
