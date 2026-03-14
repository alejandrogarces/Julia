using LinearAlgebra
using Plots
using JuMP
using Ipopt

f(y) = y^2
F(y) = -5*y
A = [1 1]
b = 1.0
c = [3;-2]
redondear(x) = round(x,digits=5)
# Modelo general
ModeloGeneral = Model(Ipopt.Optimizer)
@variable(ModeloGeneral,x[1:2].>=0)
@variable(ModeloGeneral, 10>=y >=-10)
@constraint(ModeloGeneral,A*x+[F(y)] == [b])
@objective(ModeloGeneral,Min,c'*x+f(y))
set_silent(ModeloGeneral)
set_attribute(ModeloGeneral, "print_level", 0)
optimize!(ModeloGeneral)
printstyled("Modelo completo\n", color=:red)
printstyled("-------------------------------\n",color=:red)
println(termination_status(ModeloGeneral))
println("x    = ",redondear.(value.(x)))
println("y    = ",redondear(value(y)))
println("cx   = ", redondear.(c'*value.(x)))
println("f(y) = ", redondear.(f(value(y))))



# cortes de Benders
printstyled("Cortes de Benders\n", color=:green)
printstyled("-------------------------------\n",color=:green)

function ModeloEsclavo(yc)
    model1 = Model(Ipopt.Optimizer)
    @variable(model1,x[1:2]>=0)
    @constraint(model1,restriccion,A*x+[F(yc)] == [b])
    @objective(model1,Min,c'*x)
    set_silent(model1)
    set_attribute(model1, "print_level", 0)
    optimize!(model1)
    status = termination_status(model1)        
    xk = redondear.(value.(x))
    return dual(restriccion), status, xk
end


# Modelo maestro
ModeloMaestro = Model(Ipopt.Optimizer)
set_silent(ModeloMaestro)
set_attribute(ModeloMaestro, "print_level", 0)
@variable(ModeloMaestro, 10>=y >=-10)
@variable(ModeloMaestro, θ>=-1000)
@objective(ModeloMaestro, Min, f(y)+θ)

# Cortes de Benders
println("y\t x\t\t u\t\t cx\t f(y)\t Corte")
global yk = -10
for k = 1:4   
    print([redondear.(yk)],"\t")
    local uk,status,xk = ModeloEsclavo(yk)
    if (status==LOCALLY_SOLVED)|(status==OPTIMAL)
        print(redondear.(xk),"\t")
        print(redondear.(uk),"\t")  
        print(redondear.([c'*xk]),"\t")  
        print(redondear.([f(yk)]),"\t")
        println("optimalidad")
        @constraint(ModeloMaestro, θ >= (b-F(y))'*uk[1])
        #@constraint(ModeloMaestro, θ >= f(yk))
        optimize!(ModeloMaestro)
    else
        print("Infactible\t")
        print(redondear.(uk),"\t")
        print(redondear.([c'*xk]),"\t")  
        print(redondear.([f(yk)]),"\t")
        println("factibilidad")
        @constraint(ModeloMaestro, (b-F(y))'*uk[1] <= 0)    
        optimize!(ModeloMaestro)    
    end  
    global yk = value.(y)
end
