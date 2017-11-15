include("/Users/guilhermebodin/Documents/PUC/MIPTests_Branch_Guilherme/miptests.jl")
include("/Users/guilhermebodin/Documents/PUC/MIPSolver2017_04.jl/src/branch_and_bound.jl")

function testRoutine(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Testes" begin

        solver = GurobiSolver(OutputFlag = 0)

        test1(solveMIP, solver)
        test2(solveMIP, solver)
        testSudoku(solveMIP, solver)
        test3(SolveMIP, solver)
        testInfeasibleKnapsack(solveMIP, solver)
        test_P1_Brito(solveMIP, solver)
        test_PL_Simples_Raphael(solveMIP, solver)
        test_PL_Infeasible_Brito(solveMIP, solver)
        test_PL_Unbounded_Brito(solveMIP, solver)
        test_MIP_Minimal_Brito(solveMIP, solver)
        test_MIP_Pequeno_Brito(solveMIP, solver)
        #testRobustCCUC(solveMIP, solver)
        testCaminho(solveMIP, solver)
        test3_2(solveMIP, solver)
        test3_3(solveMIP, solver)
        test_feature_selection_pequeno_viavel(solveMIP, solver)
        test_feature_selection_medio(solveMIP, solver)
        test_feature_selection_grande(solveMIP, solver)
        test_feature_selection_pequeno_inviavel(solveMIP, solver)
        teste_PL_andrew_unbounded(solveMIP, solver)
        teste_PL_andrew_viavel(solveMIP, solver)
        teste_PL_andrew_inviavel(solveMIP, solver)
        testUnboundedKnapsack(solveMIP, solver)
        testInfeasibleUC(solveMIP, solver)
        test_PL_Simples_Raphael(solveMIP, solver)
        test_PL_Infeasible_Raphael(solveMIP, solver)
        test_Minimal_UC(solveMIP, solver)
        testSudoku4x4(solveMIP, solver)

    end
end

testRoutine(SolveMIP)

#teste P1 TSP de 7 cidades
#adicionado por Guilherme Bodin
function test_P1_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste da P1 Guilherme (TSP 7 cidades)" begin
        number_of_nodes = 7

        C = [0.0    135.484   142.801   131.0     117.154   153.473   201.022
             135.484    0.0     105.546   174.003   142.425    53.0094  105.991
             142.801  105.546     0.0      87.6641   59.0      63.8905   73.5527
             131.0    174.003    87.6641    0.0      31.9061  146.932   159.201
             117.154  142.425    59.0      31.9061    0.0     115.521   132.306
             153.473   53.0094   63.8905  146.932   115.521     0.0      55.4617
             201.022  105.991    73.5527  159.201   132.306    55.4617    0.0   ]

        ans = [0.0   0.0   0.0   1.0   0.0   0.0   0.0
               1.0   0.0   0.0   0.0   0.0   0.0   0.0
               0.0   0.0   0.0   0.0   0.0   0.0   1.0
               0.0   0.0   0.0   0.0   1.0   0.0   0.0
               0.0   0.0   1.0   0.0   0.0   0.0   0.0
               0.0   1.0   0.0   0.0   0.0   0.0   0.0
               0.0   0.0   0.0   0.0   0.0   1.0   0.0]

        m = Model(solver = solver)
        @variable(m, X[i=1:number_of_nodes,j=1:number_of_nodes], Bin)
        @variable(m, u[i=:1:number_of_nodes], Int)
        for i=1:number_of_nodes
            @constraint(m,sum(X[i,j] for j=1:number_of_nodes if j!=i) == 1 )
        end
        for j=1:number_of_nodes
            @constraint(m,sum(X[i,j] for i=1:number_of_nodes if i!=j) == 1 )
        end
        for i=2:number_of_nodes
            for j=2:number_of_nodes
              if(i!=j)
                @constraint(m, u[i] - u[j] + number_of_nodes*X[i,j] <= number_of_nodes-1)
              end
            end
        end
        @objective(m, Min, sum(C[i,j]*X[i,j] for i=1:number_of_nodes, j=1:number_of_nodes))

        sol = solveMIP(m)
        @test getobjectivevalue(m) == 539.4139
        @test getvalue(X) == ans || getvalue(X) == ans'
    end
end

function test_PL_Infeasible_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste PL Infeasible Guilherme" begin
        m = Model(solver = solver)
        @variable(m, x[i=1:2])
        @constraint(m, x[1] == 6)
        @constraint(m, x[2] == 6)
        @constraint(m, x[1] + x[2] <=11)
        @objective(m, Min, x[1]+x[2])

        sol = solveMIP(m)
        @test m.ext[:status] == :Infeasible
    end
end

#teste MIP médio TSP de 20 cidades
#adicionado por Guilherme Bodin
function test_MIP_médio_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste MIP médio Guilherme (TSP 30 cidades)" begin
        number_of_nodes = 30
        srand(12)
        C = 1000*rand(30,30)

        m = Model(solver = solver)
        @variable(m, X[i=1:number_of_nodes,j=1:number_of_nodes], Bin)
        @variable(m, u[i=:1:number_of_nodes], Int)
        for i=1:number_of_nodes
            @constraint(m,sum(X[i,j] for j=1:number_of_nodes if j!=i) == 1 )
        end
        for j=1:number_of_nodes
            @constraint(m,sum(X[i,j] for i=1:number_of_nodes if i!=j) == 1 )
        end
        for i=2:number_of_nodes
            for j=2:number_of_nodes
              if(i!=j)
                @constraint(m, u[i] - u[j] + number_of_nodes*X[i,j] <= number_of_nodes-1)
              end
            end
        end
        @objective(m, Min, sum(C[i,j]*X[i,j] for i=1:number_of_nodes, j=1:number_of_nodes))

        sol = solveMIP(m)
        @test getobjectivevalue(m) ≈ 1645.8508340848819 atol = 1e-7
    end
end

#teste MIP grande TSP de 100 cidades
#adicionado por Guilherme Bodin
function test_MIP_Grande_Guilherme(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Teste MIP Grande Guilherme (TSP 100 cidades)" begin
        number_of_nodes = 100
        srand(12)
        A = 1000*rand(100,100)

        m = Model(solver = solver)
        @variable(m, X[i=1:number_of_nodes,j=1:number_of_nodes], Bin)
        @variable(m, u[i=:1:number_of_nodes], Int)
        for i=1:number_of_nodes
            @constraint(m,sum(X[i,j] for j=1:number_of_nodes if j!=i) == 1 )
        end
        for j=1:number_of_nodes
            @constraint(m,sum(X[i,j] for i=1:number_of_nodes if i!=j) == 1 )
        end
        for i=2:number_of_nodes
            for j=2:number_of_nodes
              if(i!=j)
                @constraint(m, u[i] - u[j] + number_of_nodes*X[i,j] <= number_of_nodes-1)
              end
            end
        end
        @objective(m, Min, sum(C[i,j]*X[i,j] for i=1:number_of_nodes, j=1:number_of_nodes))

        sol = solveMIP(m)
        @test getobjectivevalue(m) ≈ 1720.190204078063 atol = 1e-7
    end
end




test_P1_Guilherme(SolveMIP, GurobiSolver())
test_PL_Infeasible_Guilherme(SolveMIP, GurobiSolver())
test_MIP_médio_Guilherme(SolveMIP, GurobiSolver())
test_MIP
