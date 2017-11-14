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

function IndeciseVariable(node::Node)
    sizeX = length(node.Xrelax)
    nonInteger = Array{Float64}(sizeX)
    for i=1:sizeX
        if node.Model.colCat[i] == :Bin
            nonInteger[i] = abs(node.Xrelax[i] - 0.5)
        else
            nonInteger[i] = 2
        end
    end
    mostInDoubt = indmin(nonInteger) #get index of the most fractionary variable
    if nonInteger[mostInDoubt] == 0.5  #all variables are 0 or 1
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
    indeciseVariable = IndeciseVariable(node)
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

function UpdateBest(best::Best, node::Node, sense::Symbol, binaryVariables::Array{Int})
    if sense == :Max
        if (node.Zbound > best.Zstar && all(isinteger, node.Xrelax[binaryVariables]))
            best.Zstar = node.Zbound
            best.Xstar = node.Xrelax
        end
    else
        if (node.Zbound < best.Zstar && all(isinteger, node.Xrelax[binaryVariables]))
            best.Zstar = node.Zbound
            best.Xstar = node.Xrelax
        end
    end
end

function BinaryVariables(model::JuMP.Model)
    binaryVariableIndexes = Array{Int}(0)
    for i=1:length(model.colCat)
        if model.colCat[i] == :Bin
            push!(binaryVariableIndexes,i)
        end
    end
    return binaryVariableIndexes
end

function SolveMIP(model::JuMP.Model)
    tic()
    binaryVariableIndexes = BinaryVariables(model)
    solver = GurobiSolver(OutputFlag=0)
    iter = 0
    level = 0
    sense = getobjectivesense(model)
    best = InitializeBest(model)
    nodes = Array{Node}(1)
    nodes[1] = InitializeHeadNode(model, sense) #nó raiz
    nodes[1].Status, nodes[1].Zbound, nodes[1].Xrelax = SolveRelax(nodes[1].Model,solver)
    model.ext[:status] = nodes[1].Status
    UpdateBest(best, nodes[1], sense, binaryVariableIndexes)
    best.Visited = best.Visited + 1

    while (!isempty(nodes) && iter <= 1000)
        level = level + 1
        println(level)
        leftChild, rightChild = Branch(nodes[end])
        pop!(nodes)

        if leftChild != false
            statusLeftChild, boundLeftChild, xRelaxLeftChild = SolveRelax(leftChild, solver)
            best.Visited = best.Visited + 1
            nodeLeftChild = Node(level, leftChild, boundLeftChild, xRelaxLeftChild, statusLeftChild)
            UpdateBest(best, nodeLeftChild, sense, binaryVariableIndexes)
            if Bound(nodeLeftChild, best, sense) != false
                push!(nodes, nodeLeftChild)
            end
        end

        if rightChild != false
            statusRightChild, boundRightChild, xRelaxRightChild = SolveRelax(rightChild, solver)
            best.Visited = best.Visited + 1
            nodeRightChild = Node(level, rightChild, boundRightChild, xRelaxRightChild, statusRightChild)
            UpdateBest(best,nodeRightChild, sense, binaryVariableIndexes)
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


Cinv = 13.16
M = 200
model = Model()
@variable(model, x[i=1:2]>=0)
@variable(model, u, Bin)
@objective(model, Max, 4*x[1] + 3*x[2] - u*Cinv)
@constraint(model, 2*x[1] + 1*x[2] <= 4 +u*M)
@constraint(model, 1*x[1] + 2*x[2] <= 4 +u*M)
@constraint(model, 1*x[1] + 0.1*x[2] <= 4 +(1-u)*M)
@constraint(model, 0.4*x[1] + 1*x[2] <= 4 +(1-u)*M)

SolveMIP(model)

model.ext[:Visited]
model.ext[:status]
model.objVal
model.colVal



m = Model()
@variable(m, x[i=1:3], Bin)
@constraint(m, 6*x[1] + 5*x[2] + 5*x[3] <= 10)
@objective(m, Max, 6*x[1] + 4*x[2] + 3*x[3])

SolveMIP(m)
m.ext[:Visited]
