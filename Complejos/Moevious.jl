using Plots

function PlotTransformacion(f,puntos)
    theme(:dark)
    xmax = maximum([real.(puntos) 1])
    ymax = maximum([imag.(puntos) 1])
    xmin = minimum([real.(puntos) -1])
    ymin = minimum([imag.(puntos) -1])
    V = [xmin+ymin*1im xmax+ymin*1im xmax+ymax*1im xmin+ymax*1im]
    t = 0:0.01:1
    Vf = f.(V)
    p = plot([f.(puntos[1])],seriestype=:scatter, label="",markersize=5,aspect_ratio=:equal)
    for k = 2:length(puntos)
        plot!([f.(puntos[k])],seriestype=:scatter, label="",markersize=5,aspect_ratio=:equal)
    end
    pasox = (xmax-xmin)/10
    pasoy = (ymax-ymin)/10
    linea(a,b,λ) = a*λ + b*(1-λ) 
    for k = 0:10
        z = linea.(xmin+(ymin+k*pasoy)*1im,xmax+(ymin+k*pasoy)*1im,t)
        p = plot!(f.(z),label="",linecolor=:gray)
        z = linea.(k*pasox+xmin+ymin*1im,k*pasox+xmin+ymax*1im,t)
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
    p = annotate!(real(Vf[1]), imag(Vf[1]), string(1), :white)
    return p
end

function AnimarTransformacion(f,puntos,n)
    r(k) = 1-k^4
    anim = @animate for k ∈ -n:n
        λ = r(k/n)
        fk(z) = (1-λ)*z + λ*f(z)
        PlotTransformacion(fk,puntos)
    end
    gif(anim, "TransformacionMoevious.gif", fps = n)
end


## Ejemplos de Transformaciones de Moevious
a = 0
b = 2
c = 1
d = 2
T(z) = (a*z+b)/(c*z+d) 
puntos = [0 -d/c*0.9] .+ 0im
println("det=",a*d-b*c)
println("polo=",-d/c)
AnimarTransformacion(T,puntos,100)
