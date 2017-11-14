using JuMP, Gurobi

mutable struct Node

    Level::Int
    Model::JuMP.Model
    Z_ub::Float64
    Z_lb::Float64
    Xrelax::Array{Float64}
    Status::Symbol

end

mutable struct Best

    Zstar::Float64
    Xstar::Array{Float64}
    Visited::Int

end

function IndeciseVariable(Xrelax::Array{Float64})
    sizeX = length(Xrelax)
    nonInteger = Array{Float64}(sizeX)

    for i=1:sizeX
        nonInteger[i] = abs(Xrelax[i] - 0.5)
    end

    mostInDoubt = indmin(nonInteger) #get index of the most fractionary variable

    if nonInteger[mostInDoubt] == 0.5 #all variables are 0 or 1
        return false
    end

    return mostInDoubt
end

function Bound(node::Node, best::Best, sense::Symbol, errortolerance::Float64)

    errortolerance = 1e-4 #if 0.01% away from optimal -> current answer is optimal

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

    if abs(node.Z_ub - node.Z_lb) <= errortolerance #Bound by Optimality
        return false
    end
end

function Branch(node::Node)

    indeciseVariable = IndeciseVariable(node.Xrelax)
    if !indeciseVariable #if all variables are integer -> no branching
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
    return leftChild, rightChild
end

function InitializeHeadNode(model::JuMP.Model)
            #Level,    Model,   Z_ub,  Z_lb,           Xrelax,                        Status
    node = Node(0, copy(model), 1e10, -1e10, Array{Float64}(length(model.colUpper)), :None)
end

function InitializeBest(model::JuMP.Model)
    sense = getobjectivesense(model)
               #Zstar,              Xstar,                 Visited
    best = Best(-Inf, Array{Float64}(length(model.colUpper)),0) #initialize Best as a Max problem

    if sense == :Min
        best.Zstar = Inf
    end

    return best
end

function SolveRelax(model::JuMP.Model)


function solveMIP(model::JuMP.Model)
    InitializeBest(model)
    InitializeHeadNode(model)
    Xstar = []
    nodes = Array{node}(1)
    nodes[1] = node(0,model) #nó raiz

end

a = :Min
typeof(a)

model = Model(solver = GurobiSolver())
@variable(model, x[i=1:3], Bin)
@constraint(model, 6*x[1] + 5*x[2] + 5*x[3] <= 10)
@objective(model, Max, 6*x[1] + 4*x[2] + 3*x[3])
model
status = solve(model, relaxation = true)
X = getvalue(x)
IndeciseVariable(X)
getobjectivevalue(model)
X = model.colUpper
model.colLower
model.obj
getobjectivevalue(model)


m = copy(model)

m.colUpper

m

solveMIP(model)

model.colUpper
model.colLower

model.colCat[:]

sense = getobjectivesense(model)
typeof(sense)

nodes = Array{node}(1)
nodes[1] = node(0,model)

nodes
leftChild = node(1,model)
rightChild = node(1,model)

nodes = push!(nodes, leftChild)

v = Array{Float64}(3)

v[1]=20
v[2]=2
v[3]=4

indmin(v)

a = false

if !a
    print(1)
end

a::Float64
b = 4.5434
typeof(b)
b=Inf
a = Inf
typeof(b)

a = :Optimal

typeof(a)
