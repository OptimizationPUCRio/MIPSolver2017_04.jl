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
    Intsols::Int
    Visited::Int
end

function IndeciseVariable(node::Node)
    sizeX = length(node.Xrelax)
    nonInteger = Array{Float64}(sizeX)
    for i=1:sizeX
        if node.Model.colCat[i] == :Bin
            nonInteger[i] = abs(node.Xrelax[i] - 0.5)
        else
            nonInteger[i] = 2 #Variable is not binary
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
    if indeciseVariable == false #if all binary variables are integer -> no branching
        return false, false
    end

    leftChild = deepcopy(node.Model)
    # Forces the most indecise variable to be 0
    leftChild.colUpper[indeciseVariable] = 0
    leftChild.colLower[indeciseVariable] = 0

    rightChild = deepcopy(node.Model)
    # Forces the most indecise variable to be 0
    rightChild.colUpper[indeciseVariable] = 1
    rightChild.colLower[indeciseVariable] = 1

    return leftChild, rightChild
end

function InitializeHeadNode(model::JuMP.Model, sense::Symbol)
            #Level,    Model,   Zbound,             Xrelax,                 Status
    node = Node(0, deepcopy(model), Inf, Array{Float64}(length(model.colUpper)), :None)
    if sense == :Min
        node.Zbound = -Inf
    end
    return node
end

function InitializeBest(model::JuMP.Model)
    sense = getobjectivesense(model)
               #Zstar,         Xstar,                  Intsols    Visited
    best = Best( -Inf, Array{Float64}(length(model.colUpper)), 0 ,0) #initialize Best as a Max problem
    if sense == :Min
        best.Zstar = Inf
    end
    return best
end

function SolveRelax(model::JuMP.Model, solver::MathProgBase.AbstractMathProgSolver)
    setsolver(model, solver)
    status = solve(model, relaxation = true) #solve model relaxation
    return status, getobjectivevalue(model), model.colVal
end

function UpdateBest(best::Best, node::Node, sense::Symbol, binaryVariables::Array{Int})
    if sense == :Max
        if (node.Zbound > best.Zstar && all(isinteger, node.Xrelax[binaryVariables])) #if the bound is better and if all binary variables are integer
            best.Zstar = node.Zbound
            best.Xstar = node.Xrelax
            best.Intsols = best.Intsols + 1
        end
    else
        if (node.Zbound < best.Zstar && all(isinteger, node.Xrelax[binaryVariables]))
            best.Zstar = node.Zbound
            best.Xstar = node.Xrelax
            best.Intsols = best.Intsols + 1
        end
    end
end

function BinaryVariables(model::JuMP.Model)
    binaryVariableIndexes = Array{Int}(0)
    for i=1:length(model.colCat)
        if model.colCat[i] == :Bin
            push!(binaryVariableIndexes,i) #gets index of all binary variables
        end
    end
    return binaryVariableIndexes
end

function SolveMIP(model::JuMP.Model)
    tic()
    iter = 0
    maxiter = 1000
    level = 0
    solver = GurobiSolver(OutputFlag=0)
    sense = getobjectivesense(model)
    best = InitializeBest(model)
    nodes = Array{Node}(1)
    binaryVariableIndexes = BinaryVariables(model) #gets index of all binary variables
    nodes[1] = InitializeHeadNode(model, sense) #root node
    nodes[1].Status, nodes[1].Zbound, nodes[1].Xrelax = SolveRelax(nodes[1].Model,solver)
    model.ext[:status] = nodes[1].Status
    UpdateBest(best, nodes[1], sense, binaryVariableIndexes) #updates best answer
    best.Visited = best.Visited + 1

    while (!isempty(nodes) && iter <= maxiter)
        level = level + 1
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

    if (iter >= maxiter && all(isinteger, best.Xstar[binaryVariableIndexes]))
        model.ext[:status] = :SubOptimal
    elseif (iter >= maxiter && ~(all(isinteger, best.Xstar[binaryVariableIndexes])))
        mode.ext[:status] = :NotResolved
    end

    model.ext[:nodes] = best.Visited
    model.ext[:intsols] = best.Intsols
    model.colVal = best.Xstar
    model.objVal = best.Zstar
    model.objBound = best.Zstar
    time = toc()
    model.ext[:time] = time
end
