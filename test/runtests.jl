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
        testMinimalTSP(solveMIP)
    end
end

testRoutine(SolveMIP)
