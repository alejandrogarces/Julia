using Agents
using Plots
using Random
# 1. Definición del Agente
@agent struct Car(ContinuousAgent{2, Float64})
    is_aggressive::Bool
    speed::Float64
end
# 2. Inicialización del Modelo
function initialize_model(; N=80, seed=42)
    my_rng = MersenneTwister(seed)    
    space = ContinuousSpace((500.0, 500.0); periodic = false)    
    model = StandardABM(Car,space;agent_step! = car_step!,rng = my_rng )
    for i in 1:N
        is_aggr = rand(abmrng(model)) > 0.8
        x_start = is_aggr ? rand(abmrng(model), 0.0:8.0) : rand(abmrng(model), 12.0:14.0)
        y_start = rand(abmrng(model), 100.0:500.0)
        pos = (x_start, y_start)
        vel = (0.0, -1.5)
        add_agent!(pos, model, vel, is_aggr, 1.5)
    end
    return model
end
function car_step!(car, model)
    x, y = car.pos
    target_speed = car.speed
    rng = abmrng(model)
    noise_x = rand(rng, -0.08:0.01:0.08)
    noise_y = rand(rng, -0.05:0.01:0.05)
        min_dist_ahead = Inf
    
    for other in nearby_agents(car, model, 2.0)
        # El otro auto está adelante si su coordenada Y es menor (bajamos hacia 0)
        if other.pos[2] < y
            # Distancia euclidiana combinada en X e Y
            dist = hypot(other.pos[1] - x, y - other.pos[2])
            if dist < 2.0 && dist < min_dist_ahead
                min_dist_ahead = dist
            end
        end
    end
    
    # Si hay un vehículo dentro del radio estricto de 2 unidades, se detiene por completo
    if min_dist_ahead < 2.0
        target_speed = 0.0
    end
    
    # Cálculo de dirección base hacia abajo con ruido estocástico
    dx = 0.0 + noise_x
    dy = -target_speed + noise_y
    
    # Zona de inflexión / embudo (y < 40): giro de 30° a la derecha y reducción de carriles
    if y < 40.0
        dx = target_speed * sin(deg2rad(60)) + noise_x
        dy = -target_speed * cos(deg2rad(60)) + noise_y
        
        # Los agresivos intentan colarse forzando el paso a la derecha en X
        if car.is_aggressive && x < 50.0
            dx += 1.0
        end
    elseif car.is_aggressive && y < 60.0
        if target_speed > 0.0
            dy = -target_speed * 1.3 + noise_y
        end
    end

    # Actualizar velocidad y mover al agente
    car.vel = (dx, dy)
    move_agent!(car, model, 1.0)
end
# 4. Creación y Ejecución de la Simulación
model = initialize_model(N=800)
theme(:dark)
# 5. Generación del Video / GIF usando Plots.jl
println("Generando simulación...")
anim = @animate for t in 1:300
    step!(model, 1)
    # Extraer posiciones y tipos de conductores
    xs = [a.pos[1] for a in allagents(model)]
    ys = [a.pos[2] for a in allagents(model)]
    colors = [a.is_aggressive ? :red : :blue for a in allagents(model)]
    
    # Graficar límites de la carretera
    plot([15, 15, 50], [100, 40, 28], color=:gray, label="", lw=2,axis=false, ticks=false)
    plot!([ 0,  0, 50], [100, 15, 15], color=:gray, label="", lw=2,axis=false, ticks=false)

    # Graficar vehículos (Azul = Precavido, Rojo = Agresivo)
    scatter!(xs, ys, color=colors, ms=2, legend=false,
             xlims=(0, 50), ylims=(0, 100), 
             title="Entrega de medallas MMC",aspect_ratio=:equal)
end

# Guardar animación final
println("¡Animación guardada como 'simulacion.gif'!")
gif(anim, "simulacion.gif", fps = 15)

# Crear un video
using Pkg
Pkg.add("FFMPEG")
mp4(anim, "simulacion_medallas_mmc.mp4", fps = 15)
