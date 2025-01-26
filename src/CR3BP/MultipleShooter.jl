"""
CR3BP multiple shooter wrapper

Author: Jonathan Richmond
C: 9/8/22
U: 1/22/25
"""

import LinearAlgebra, StaticArrays
import MBD: CR3BPMultipleShooter

export solveUpdateEquation

"""
    solve!(multipleShooter, initialGuess)

Return converged solution

# Arguments
- `multipleShooter::CR3BPMultipleShooter`: CR3BP multiple shooter object
- `initialGuess::CR3BPMultipleShooterProblem`: CR3BP unsolved multiple shooter problem
"""
function solve!(multipleShooter::CR3BPMultipleShooter, initialGuess::MBD.CR3BPMultipleShooterProblem)
    buildProblem!(initialGuess)
    solutionInProgress::MBD.CR3BPMultipleShooterProblem = deepClone(initialGuess)
    numFreeVariables::Int16 = Int16(getNumFreeVariables!(solutionInProgress))
    ((numFreeVariables == Int16(0)) || (getNumConstraints(solutionInProgress) == 0)) ? (return solutionInProgress) : multipleShooter.recentIterationCount = 0
    while !isConverged(multipleShooter.convergenceCheck, solutionInProgress) && (multipleShooter.recentIterationCount < multipleShooter.maxIterations)
        if multipleShooter.recentIterationCount > 0
            freeVariableStep::StaticArrays.SVector{Int64(numFreeVariables), Float64} = solveUpdateEquation(multipleShooter, solutionInProgress)
            freeVariableVector::Vector{Float64} = getFreeVariableVector!(solutionInProgress)+freeVariableStep
            setFreeVariableVector!(solutionInProgress, freeVariableVector)
        end
        multipleShooter.printProgress && println("Iteration $(multipleShooter.recentIterationCount): ||F|| = $(LinearAlgebra.norm(getConstraintVector!(solutionInProgress)))")
        multipleShooter.recentIterationCount += 1
    end
    
    isConverged(multipleShooter.convergenceCheck, solutionInProgress) ? (return solutionInProgress) : throw(ErrorException("Corrections algorithm could not converge"))
end

"""
    solveUpdateEquation(multipleShooter, multipleShooterProblem)

Return free variable vector update

# Arguments
- `multipleShooter::CR3BPMultipleShooter`: CR3BP multiple shooter object
- `multipleShooterProblem::CR3BPMultipleShooterProblem`: CR3BP multiple shooter problem object
"""
function solveUpdateEquation(multipleShooter::CR3BPMultipleShooter, multipleShooterProblem::MBD.CR3BPMultipleShooterProblem)
    for generator::MBD.AbstractUpdateGenerator in multipleShooter.updateGenerators
        canGenerateUpdate(generator, multipleShooterProblem) && (return getFullUpdate(generator, multipleShooterProblem))
    end
    throw(ErrorException("No update generators were successful"))
end
