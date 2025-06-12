"""
BCR4BP P1-P2 multiple shooter wrapper

Author: Jonathan Richmond
C: 4/9/25
"""

import LinearAlgebra, StaticArrays
import MBD: BCR4BP12MultipleShooter

export solveUpdateEquation

"""
    solve!(multipleShooter, initialGuess)

Return converged solution

# Arguments
- `multipleShooter::BCR4BP12MultipleShooter`: BCR4BP P1-P2 multiple shooter object
- `initialGuess::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 unsolved multiple shooter problem
"""
function solve!(multipleShooter::BCR4BP12MultipleShooter, initialGuess::MBD.BCR4BP12MultipleShooterProblem)
    buildProblem!(initialGuess)
    solutionInProgress::MBD.BCR4BP12MultipleShooterProblem = deepClone(initialGuess)
    numFreeVariables::Int16 = Int16(getNumFreeVariables!(solutionInProgress))
    ((numFreeVariables == Int16(0)) || (getNumConstraints(solutionInProgress) == 0)) ? (return solutionInProgress) : multipleShooter.recentIterationCount = 0
    while !isConverged(multipleShooter.convergenceCheck, solutionInProgress) && (multipleShooter.recentIterationCount < multipleShooter.maxIterations)
        if multipleShooter.recentIterationCount > 0
            freeVariableStep::StaticArrays.SVector{Int64(numFreeVariables), Float64} = solveUpdateEquation(multipleShooter, solutionInProgress)
            println(freVariableStep)
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
- `multipleShooter::BCR4BP12MultipleShooter`: BCR4BP P1-P2 multiple shooter object
- `multipleShooterProblem::BCR4BP12MultipleShooterProblem`: BCR4BP P1-P2 multiple shooter problem object
"""
function solveUpdateEquation(multipleShooter::BCR4BP12MultipleShooter, multipleShooterProblem::MBD.BCR4BP12MultipleShooterProblem)
    for generator::MBD.AbstractUpdateGenerator in multipleShooter.updateGenerators
        canGenerateUpdate(generator, multipleShooterProblem) && (return getFullUpdate(generator, multipleShooterProblem))
    end
    throw(ErrorException("No update generators were successful"))
end
