############################################
####  Mapas varios en variable compleja ####
####  Ejemplos simples que demuestran   ####
####  como se transforma el espacio     ####
####  con diferentes mapas complejos    ####
####  Por:  Alejandro Garces            ####
####        19/12/2025                  ####
############################################

using Plots

function PlotTransformacion(V,f)
    theme(:dark)
    t = 0:0.01:1
    Vf = f.(V)
    p = plot([f(0+0im)+0im],seriestype=:scatter, label="",markersize=5)
    xmax = maximum(real.(V))
    ymax = maximum(imag.(V))
    xmin = minimum(real.(V))
    ymin = minimum(imag.(V))
    dx = (xmax-xmin)/10
    dy = (ymax-ymin)/10
    for k = 0:10
        z = ((ymax-ymin)*t.+ymin)*1im.+(xmin+k*dx)
        p = plot!(f.(z),label="",linecolor=:gray)
        z = ((xmax-xmin)*t.+xmin).+(ymin+k*dy)*1im
        p = plot!(f.(z),label="",linecolor=:gray)
    end
    n = length(V)
    for k = 2:n            
        z = (V[k-1]-V[k])*t.+V[k]        
        w = f.(z)
        p=plot!(w,label="",linecolor=:orange, linewidth=3)
        p = annotate!(real(Vf[k]), imag(Vf[k]), string(k), :white)
    end
    z = (V[n]-V[1])*t.+V[1]
    w = f.(z)
    p=plot!(w,label="",linecolor=:orange, linewidth=3)
    p = plot!(Vf[:],seriestype=:scatter, label="", markersize=12)
    p = annotate!(real(Vf[1]), imag(Vf[1]), string(1), :white)
    return p
end

function AnimarTransformacion(V,f,n)
    anim = @animate for k ∈ 0:n
        λ = abs(2*k/n-1)
        fk(z) = λ*z + (1-λ)*f(z)
        PlotTransformacion(V,fk)
    end
    gif(anim, "TRANSFORMACION.gif", fps = n)
end

## Ejemplos de puntos clave
#V = randn(4)+randn(4)*1im
#V = [1+0.5im 1.5+0.5im 1.5+1im 1+1im]
#V = [-0.5-0.5im 0.5-0.5im 0.5+0.5im -0.5+0.5im]
V = [-0.5+0.0im 0.0-0.5im 0.5+0.0im 0.0+0.5im]

## Ejemplos de mapas en variable compleja
#f(z) = exp(z)
#f(z) = exp(-z)
#f(z) = z*z
f(z) = z^3
AnimarTransformacion(V,f,200)
