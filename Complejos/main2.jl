############################################
####  Mapas varios en variable compleja ####
####  Ejemplos simples que demuestran   ####
####  como se transforma el espacio     ####
####  con diferentes mapas complejos    ####
####  Por:  Alejandro Garces            ####
####        19/12/2025                  ####
############################################
using Plots

# todos los graficos estaran en un [-1,1]x[-1,1]
function PlotTransformacion(f,ξ)
    theme(:dark)
    V = [-1-1im 1-1im 1+1im -1+1im]
    t = 0:0.01:1
    Vf = f.(V)
    p = plot([f(0+0im)+0im],seriestype=:scatter, label="",markersize=5,aspect_ratio=:equal)
    for k = 0:10
        z = (2*t.-1)*1im.+(-1+k*0.2)
        p = plot!(f.(z),label="",linecolor=:gray)
        z = (2*t.-1).+(-1+k*0.2)*1im
        p = plot!(f.(z),label="",linecolor=:gray)
    end
    for k = 2:4          
        z = (V[k-1]-V[k])*t.+V[k]        
        w = f.(z)
        p=plot!(w,label="",linecolor=:orange, linewidth=2)
        p = annotate!(real(Vf[k]), imag(Vf[k]), string(k), :white)
    end
    z = (V[4]-V[1])*t.+V[1]
    w = f.(z)
    p=plot!(w,label="",linecolor=:orange, linewidth=2)
    p = plot!(Vf[:],seriestype=:scatter, label="", markersize=12)
    p = plot!(f.(0.9*(1.0.-ξ*t).*exp.(2*pi*t*1im)), linewidth=3, color=:white,label="")
    p = annotate!(real(Vf[1]), imag(Vf[1]), string(1), :white)
    return p
end

function AnimarTransformacion(f,ξ,n)
    anim = @animate for k ∈ 0:n
        λ = abs(2*k/n-1)
        fk(z) = λ*z + (1-λ)*f(z)
        PlotTransformacion(fk,ξ)
    end
    gif(anim, "TransformacionTrayectoria.gif", fps = n)
end


## Ejemplos de mapas en variable compleja
f(z) = exp(z)
f(z) = z^3
AnimarTransformacion(f,1,200)
