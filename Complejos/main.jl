using LinearAlgebra
using Plots

function PlotBase(V)
    theme(:dark)
    xmax = maximum(abs.(real.(V)))
    ymax = maximum(abs.(imag.(V)))
    dx = xmax/5
    dy = ymax/5 
    p = plot([0],[0],seriestype=:scatter, label="",markersize=5)
    for k = 0:10        
        p = plot!([-xmax+k*dx;-xmax+k*dx],[-ymax;ymax],label="",linecolor=:gray)
        p = plot!([-xmax;xmax],[-ymax+k*dy;-ymax+k*dy],label="",linecolor=:gray)
    end
    p = plot!(V[:],seriestype=:scatter, label="",markersize=12)
    p = annotate!(real(V[1]), imag(V[1]), string(1), :white)
    t = 0:0.01:1
    n = length(V)
    for k = 2:n    
        z = (V[k-1]-V[k])*t.+V[k]
        p=plot!(z,label="",linecolor=:orange)
        p = annotate!(real(V[k]), imag(V[k]), string(k), :white)
    end
    z = (V[n]-V[1])*t.+V[1]
    p=plot!(z,label="",linecolor=:orange)
    return p
end

function PlotTransformacion(V,f)
    theme(:dark)
    t = 0:0.01:1
    Vf = f.(V)
    p = plot([f(0+0im)+0im],seriestype=:scatter, label="",markersize=5)
    xmax = maximum(abs.(real.(V)))
    ymax = maximum(abs.(imag.(V)))
    dx = xmax/5
    dy = ymax/5 
    for k = 0:10
        z = (-2*ymax*t.+ymax)*1im.+(-xmax+k*dx)
        p = plot!(f.(z),label="",linecolor=:gray)
        z = (-2*xmax*t.+xmax).+(-ymax+k*dy)*1im
        p = plot!(f.(z),label="",linecolor=:gray)
    end


    p = plot!(Vf[:],seriestype=:scatter, label="", markersize=12)
    p = annotate!(real(Vf[1]), imag(Vf[1]), string(1), :white)
    
    n = length(V)
    for k = 2:n            
        z = (V[k-1]-V[k])*t.+V[k]        
        w = f.(z)
        p=plot!(w,label="",linecolor=:orange)
        p = annotate!(real(Vf[k]), imag(Vf[k]), string(k), :white)
    end
    z = (V[n]-V[1])*t.+V[1]
    w = f.(z)
    p=plot!(w,label="",linecolor=:orange)
    return p
end

V = [1+0.5im 1.5+0.5im 1.5+1im 1+1im]
#f(z) = exp(z)
f(z) = z*z
p = PlotBase(V)
q = PlotTransformacion(V,f)
pq=plot(p,q,layout=(1,2),size=(800,400))
display(pq)
savefig("plot2.png")
