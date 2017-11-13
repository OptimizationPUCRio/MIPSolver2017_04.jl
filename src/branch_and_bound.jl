using JuMP, Gurobi

type Node
    Level::Int
    Model::JuMP.Model
    Zbound::Float64
    Xrelax::Array{Float64}
    Status::Symbol
end

type Best
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

function Bound(node::Node, best::Best, sense::Symbol)
    errortolerance = 1e-4 #if 0.01% away from optimal -> current answer is optimal
    if node.Status != :Optimal #Bound by Infeasabilty
        return false
    end
    if sense == :Max #Bound by Limit
        if node.Zbound < best.Zstar
            return false
        end
    else
        if node.Zbound > best.Zstar
            return false
        end
    end
    if abs(node.Zbound - best.Zstar) <= errortolerance #Bound by Optimality
        return false
    end
end

function Branch(node::Node)
    indeciseVariable = IndeciseVariable(node.Xrelax)
    if indeciseVariable == false #if all variables are integer -> no branching
        return false, false
    end

    leftChild = copy(node.Model)
    # Forces the most indecise variable to be 0
    leftChild.colUpper[indeciseVariable] = 0
    leftChild.colLower[indeciseVariable] = 0

    rightChild = copy(node.Model)
    # Forces the most indecise variable to be 0
    rightChild.colUpper[indeciseVariable] = 1
    rightChild.colLower[indeciseVariable] = 1

    return leftChild, rightChild
end

function InitializeHeadNode(model::JuMP.Model, sense::Symbol)
            #Level,    Model,   Zbound,             Xrelax,                 Status
    node = Node(0, copy(model), Inf, Array{Float64}(length(model.colUpper)), :None)
    if sense == :Min
        node.Zbound = -Inf
    end
    return node
end

function InitializeBest(model::JuMP.Model)
    sense = getobjectivesense(model)
               #Zstar,              Xstar,                 Visited
    best = Best( -Inf, Array{Float64}(length(model.colUpper)),0) #initialize Best as a Max problem
    if sense == :Min
        best.Zstar = Inf
    end
    return best
end

function SolveRelax(model::JuMP.Model, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    setsolver(model, solver)
    status = solve(model, relaxation = true)
    return status, getobjectivevalue(model), model.colVal
end

function UpdateBest(best::Best, node::Node, sense::Symbol)
    if sense == :Max
        if (node.Zbound > best.Zstar && isinteger(node.Xrelax))
            best.Zstar = node.Zbound
            best.Xstar = node.Xrelax
        end
    else
        if (node.Zbound < best.Zstar && isinteger(node.Xrelax))
            best.Zstar = node.Zbound
            best.Xstar = node.Xrelax
        end
    end
end

function SolveMIP(model::JuMP.Model, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    tic()
    iter = 0
    level = 0
    sense = getobjectivesense(model)
    best = InitializeBest(model)
    nodes = Array{Node}(1)
    nodes[1] = InitializeHeadNode(model, sense) #nó raiz
    nodes[1].Status, nodes[1].Zbound, nodes[1].Xrelax = SolveRelax(nodes[1].Model,solver)
    UpdateBest(best,nodes[1],sense)
    best.Visited = best.Visited + 1

    while (!isempty(nodes) && iter <= 1000)
        level = level + 1
        leftChild, rightChild = Branch(nodes[end])
        pop!(nodes)

        if leftChild != false
            statusLeftChild, boundLeftChild, xRelaxLeftChild = SolveRelax(leftChild, solver)
            best.Visited = best.Visited + 1
            nodeLeftChild = Node(level, leftChild, boundLeftChild, xRelaxLeftChild, statusLeftChild)
            UpdateBest(best, nodeLeftChild, sense)
            if Bound(nodeLeftChild, best, sense) != false
                push!(nodes, nodeLeftChild)
            end
        end

        if rightChild != false
            statusRightChild, boundRightChild, xRelaxRightChild = SolveRelax(rightChild, solver)
            best.Visited = best.Visited + 1
            nodeRightChild = Node(level, rightChild, boundRightChild, xRelaxRightChild, statusRightChild)
            UpdateBest(best,nodeRightChild, sense)
            if Bound(nodeRightChild, best, sense) != false
                push!(nodes, nodeRightChild)
            end
        end
        iter = iter + 1
    end
    time = toc()

    m.ext[:Visited] = best.Visited
    m.ext[:Time] = time
    model.colVal = best.Xstar
    model.objVal = best.Zstar
end

SolveMIP(m, GurobiSolver())
m.objVal
m.colVal

m = Model()
@variable(m, x[i=1:3], Bin)
@constraint(m, 6*x[1] + 5*x[2] + 5*x[3] <= 10)
@objective(m, Max, 6*x[1] + 4*x[2] + 3*x[3])


n = 9
model = Model()
@variable(model, x[i in 1:n, j in 1:n, k in 1:n], Bin)

fixas = [(1,3,4), (1,5,6), (1,9,2), (2,1,8), (2,3,5), (2,6,2), (2,8,3),
        (3,5,3), (3,8,6), (4,2,2), (4,3,8), (5,6,4), (6,1,7), (6,5,5),
        (6,9,9), (7,3,2), (7,6,1), (8,2,7), (8,5,4), (8,7,9), (8,8,5),
        (9,1,6), (9,8,4)]
for idx in fixas
    @constraint(model, x[idx...] == 1)
end
@constraint(model, [j in 1:n, k in 1:n], sum(x[:,j,k]) == 1)
@constraint(model, [i in 1:n, k in 1:n], sum(x[i,:,k]) == 1)
@constraint(model, [i in 1:n, j in 1:n], sum(x[i,j,:]) == 1)
@constraint(model, [p in [0,3,6], q in [0,3,6], k in 1:n], sum(sum(x[i+p,j+q,k] for i in 1:3) for j in 1:3) == 1)
@objective(model, Min, 0)

SolveMIP(model, GurobiSolver())

sum(model.colVal)


model.colVal

model = Model()
@variable(model, x[i=1:2]>=0)
@variable(model, u, Bin)
@objective(model, Max, 4*x[1] + 3*x[2] - u*Cinv)

@constraint(model, 2*x[1] + 1*x[2] <= 4 +u*M)
@constraint(model, 1*x[1] + 2*x[2] <= 4 +u*M)

@constraint(model, 1*x[1] + 0.1*x[2] <= 4 +(1-u)*M)
@constraint(model, 0.4*x[1] + 1*x[2] <= 4 +(1-u)*M)

SolveMIP(model, GurobiSolver())
