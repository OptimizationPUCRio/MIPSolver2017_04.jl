include("/Users/guilhermebodin/Documents/PUC/MIPTests_Branch_Guilherme/miptests.jl")
include("/Users/guilhermebodin/Documents/PUC/MIPSolver2017_04.jl/src/branch_and_bound.jl")

function testRoutine(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Testes" begin
        test1(solveMIP)
        test2(solveMIP)
        testSudoku(solveMIP)
        testSudoku4x4(solveMIP)
        testInfeasibleKnapsack(solveMIP)
        testInfeasibleUC(solveMIP)
        testUnboundedKnapsack(solveMIP)
        test_MIP_Minimal_Brito(solveMIP)
        test_PL_Unbounded_Brito(solveMIP)
        test_PL_Infeasible_Brito(solveMIP)
        test_PL_Simples_Raphael(solveMIP)
        test_Minimal_UC(solveMIP)
        test_PL_Infeasible_Raphael(solveMIP)
        teste_PL_andrew_inviavel(solveMIP)
        test3_2(solveMIP)
        test3_3(solveMIP)
        testMinimalTSP(solveMIP)
    end
end

testRoutine(SolveMIP)


#teste TSP 4 cidades
#Adicionado por Guilherme Bodin
function testMinimalTSP(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste TSP 4 cidades" begin
        AdjacencyMatrix = [0 1 2 3
                           1 0 2 3
                           2 2 0 3
                           3 3 3 0]

        number_of_nodes = 4

        del = [1 2 3
               1 4 5
               2 4 6
               3 5 6]

        S = [1 0 1 0 1 0 1 0 1 0 1 0 1 0
             0 1 1 0 0 1 1 0 0 1 1 0 0 1
             0 0 0 1 1 1 1 0 0 0 0 1 1 1
             0 0 0 0 0 0 0 1 1 1 1 1 1 1]

        C = [1 2 3 2 3 3]'

         function edges_no_Subconjunto(pos,Subsets,number_of_nodes) #Seleciona todas as arestas possíveis em uma dada partição do conjunto de vértices
           A = zeros(1,2)
           for i=1:number_of_nodes
             if Subsets[i,pos]==1
               for j=1:number_of_nodes
                 if Subsets[j,pos]==0
                   A = vcat(A,[i j])
                 end
               end
             end
           end
           return A
         end

         function ida_e_volta(A) # Recebe uma matriz A n x 2 e devolve uma matriz result 2*n x 2 com os elementos de A nas linhas 1 até n e o espelhamento dos elementos de A de n+1 até 2*n
           A_aux = hcat(A,A[:,1])
           A_aux = A_aux[:,2:3]
           result = vcat(A,A_aux)
           return result
         end

         function selecao_edges(e,edge_do_ciclo)
           E = []
           for i=1:size(e,1)
             for j=1:size(edge_do_ciclo,1)
               if edge_do_ciclo[j,1] == e[i,1]
                 if edge_do_ciclo[j,2] == e[i,2]
                   E = vcat(E,i)
                 end
               end
             end
           end
           return E
         end
        model = Model(solver = solver)
        @variable(model, Y[i=1:number_of_edges],Bin)
        for i=1:number_of_nodes
          @constraint(model, sum(Y[del[i,j]] for j=1:number_of_nodes-1) == 2) # Garante que todas as cidades conectadas a 2 estradas (uma de entrada e outra de saída) no caminho ótimo
        end
        for i=1:size(S,2)
          edge_do_ciclo = ida_e_volta(edges_no_Subconjunto(i,S,number_of_nodes))
          E = selecao_edges(e,edge_do_ciclo)
          @constraint(model, sum(Y[E[j]] for j=1:size(E,1)) >= 2)
        end
        @objective(model,Min,sum(C[i]*Y[i] for i=1:number_of_edges))

        sol = solveMIP(model)

        @test model.ext[:status] == :Optimal
        @test model.colVal == [1, 1, 0, 0, 1, 1]
        @test model.objVal == 9.0
    end
end
