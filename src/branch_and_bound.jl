using JuMP, Gurobi

mutable struct Node

    Level::Int
    Model::JuMP.Model
    Z_ub::Float64
    Z_lb::Float64
    Xrelax::Vector{Float64}
    Status

end

mutable struct Best

    Zstar::Float64
    Xstar::Vector{Float64}

end

function IndeciseVariable(Xrelax::Vector{Float64})
    sizeX = lenght(Xrelax)
    nonInteger = Vector{Float64}(sizeX)

    for i=1:sizeX
        nonInteger[i] = abs(node.Xrelax[i] - 0.5)
    end

    mostInDoubt = indmin(nonInteger)

    if nonInteger[mostInDoubt] == 0.5 #all variables are 0 or 1
        return false
    end

    return mostInDoubt
end

function Bound(node::Node, best::Best, sense::Symbol)

    if node.status != :Optimal #Bound by Infeasabilty
        return false
    end

    if sense == :Max #Bound by Limit
        if node.Z_lb > best.Zstar
            return false
        end
    else
        if node.Z_ub < best.Zstar
            return false
        end
    end

    if node.Z_ub == node.Zlb #Bound by Optimality
        return false
    end
end

function Branch(node::Node)

    indeciseVariable = IndeciseVariable(node.Xrelax)
    if !indeciseVariable #if all variables are integer no branching
        return false
    end

    leftChild = copy(node.Model)
    # Forces the most indecise variable to be 0
    leftChild.colUpper[indeciseVariable] == 0
    leftChild.colLower[indeciseVariable] == 0

    rightChild = copy(node.Model)
    # Forces the most indecise variable to be 0
    rightChild.colUpper[indeciseVariable] == 1
    rightChild.colLower[indeciseVariable] == 1

    node.Level = node.Level + 1

end

function solveMIP(model::JuMP.Model)
    sense = getobjectivesense(model)
    if sense == :Min
        Zstar = Inf
    else
        Zstar = -Inf
    end
    Xstar = []
    nodes = Vector{node}(1)
    nodes[1] = node(0,model) #nó raiz

end

a = :Min
typeof(a)

model = Model()
@variable(model,x[i=1:2])
@constraint(model, x[1]+0.4*x[2] <= 4)
@constraint(model, 0.1x[1]+x[2] <= 4)
@objective(model, Max, 4*x[1]+3*x[2])
model
m = copy(model)

m

solveMIP(model)

model.colUpper
model.colLower

model.colCat[:]

sense = getobjectivesense(model)
typeof(sense)

nodes = Vector{node}(1)
nodes[1] = node(0,model)

nodes
leftChild = node(1,model)
rightChild = node(1,model)

nodes = push!(nodes, leftChild)

v = Vector{Float64}(3)

v[1]=20
v[2]=2
v[3]=4

indmin(v)

a = false

if !a
    print(1)
end
