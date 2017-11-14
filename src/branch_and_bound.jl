using JuMP, Gurobi

mutable struct Node
    Level::Int
    Model::JuMP.Model
    Zbound::Float64
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

function SolveRelax(model::JuMP.Model, solver::MathProgBase.AbstractMathProgSolver)
    setsolver(model, solver)
    status = solve(model, relaxation = true)
    return status, getobjectivevalue(model), model.colVal
end

function UpdateBest(best::Best, node::Node, sense::Symbol)
    if sense == :Max
        if (node.Zbound > best.Zstar && all(isinteger, node.Xrelax))
            best.Zstar = node.Zbound
            best.Xstar = node.Xrelax
        end
    else
        if (node.Zbound < best.Zstar && all(isinteger, node.Xrelax))
            best.Zstar = node.Zbound
            best.Xstar = node.Xrelax
        end
    end
end


function SolveMIP(model::JuMP.Model)
    tic()
    solver = GurobiSolver(OutputFlag=0)
    iter = 0
    level = 0
    sense = getobjectivesense(model)
    best = InitializeBest(model)
    nodes = Array{Node}(1)
    nodes[1] = InitializeHeadNode(model, sense) #nó raiz
    nodes[1].Status, nodes[1].Zbound, nodes[1].Xrelax = SolveRelax(nodes[1].Model,solver)
    model.ext[:status] = nodes[1].Status
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

    model.ext[:Visited] = best.Visited
    model.ext[:Time] = time
    model.colVal = best.Xstar
    model.objVal = best.Zstar
end
