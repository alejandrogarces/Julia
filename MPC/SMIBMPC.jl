# MPC para un sistema idealizado de dos barras
# alejandro.garces@utp.edu.co

## Instalar paquekes
paquetes = ["Plots", "LinearAlgebra", "JuMP", "Ipopt"]
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


## Parametros
xL = 0.01    # impedancia equivalente
Hm = 3	     # inercia de la maquina
Δt = 1/60/2  # paso
np = 50      # numero de pasos a simular
wB = 2*pi*60 # frecuencia nominal	
pref = 5     # potencia de referencia
nh = 1       # horizonte

## Modelo
function SMIB(ω,δ,pm)	
	dω = 1/(2*Hm)*(pm-1/xL*sin(δ))
	dδ = wB*(ω-1)
	return dω,dδ
end

## Modelo de optimizacion usando punto medio

q1 = 1    # Factor de penalizacion en w
q2 = 1    # Factor de penalizacion en delta
r1 = 0.01 # Factor de penalizacion en p
qf = 1E6 # Factor de penalizacion final 
aref = asin(pref*xL)
model = Model(Ipopt.Optimizer)
@variable(model, xw[1:(nh+1)]) # frecuencia
@variable(model, xa[1:(nh+1)]) # angulo
@variable(model, up[1:nh]) # potencia activa
@objective(model, Min, q1*sum((xw[1:nh].-1).^2)+q2*sum((xa[1:nh].-aref).^2)+r1*sum((up.-pref).^2) + qf*q1*(xw[nh+1]-1)^2 + qf*q2*(xa[nh+1]-aref)^2)
for k = 1:nh
	xmw = (xw[k+1]+xw[k])/2 # punto medio
 	xma = (xa[k+1]+xa[k])/2 # punto medio
	dwm,dam = SMIB(xmw,xma,up[k])
	@constraint(model,xw[k+1] == xw[k] + Δt*dwm)
    @constraint(model,xa[k+1] == xa[k] + Δt*dam)
    @constraint(model, up[k] <= 10)
	@constraint(model, xw[k+1] <= 60.10/60)
	@constraint(model, xw[k+1] >= 59.90/60)
end
@constraint(model,w0, xw[1] == 0)  # condiciones iniciales
@constraint(model,a0, xa[1] == 0) # condiciones iniciales
set_silent(model)
#println(model)

## Control
function u_control(ωk,δk)#
	set_normalized_rhs(w0, ωk) 
	set_normalized_rhs(a0, δk) 
	optimize!(model)			
	return value(up[1])
end





## Simulacion usando RungeKutta
function simulacion() 
	ω = ones(np)
	δ = zeros(np)
	p = ones(np)*pref
	for k = 1:np-1
		p[k] = u_control(ω[k],δ[k])		
		fω1,fδ1 = SMIB(ω[k],δ[k],p[k])
		fω2,fδ2 = SMIB(ω[k]+Δt/2*fω1,δ[k]+Δt/2*fδ1,p[k])
		fω3,fδ3 = SMIB(ω[k]+Δt/2*fω2,δ[k]+Δt/2*fδ2,p[k])
		fω4,fδ4 = SMIB(ω[k]+Δt*fω3,δ[k]+Δt*fδ3,p[k])
		ω[k+1] = ω[k] + Δt/6*(fω1+2*fω2+2*fω3+fω4)
		δ[k+1] = δ[k] + Δt/6*(fδ1+2*fδ2+2*fδ3+fδ4)
	end
	p1 = plot(ω*60,label="Frecuencia",linecolor=:lawngreen,linewidth=2)
	p2 = plot(δ,label="Angulo",linecolor=:orange,linewidth=2)
	p3 = plot(p,label="Potencia", linetype=:steppre,linecolor=:steelblue1,linewidth=2)
	plt = plot(p1,p2,p3,layout=(3,1),size=(1200,1000))
	return plt
end

plt = simulacion()
