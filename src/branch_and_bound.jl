using JuMP

mutable struct node

    nivel::Int
    model::JuMP.Model

end

function solveMIP(model::JuMP.Model)

    if getobjectivesense(model) == :Min
        zstar = Inf
    else
        zstar = -Inf
    end
    xstar = []
    nodes = Vector{node}(1)
    nodes[1] = node(0,model) #nó raiz




end




model = Model()
@variable(model,x[i=1:2])
@constraint(model, x[1]+0.4*x[2] <= 4)
@constraint(model, 0.1x[1]+x[2] <= 4)
@objective(model, Max, 4*x[1]+3*x[2])
model
solveMIP(model)

model.colUpper
model.colLower

model.colCat[:]

getobjectivesense(model)


nodes = Vector{node}(1)
nodes[1] = node(0,model)

nodes
leftChild = node(1,model)
rightChild = node(1,model)

nodes = push!(nodes, leftChild)
