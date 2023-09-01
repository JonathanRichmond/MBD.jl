"""
Multiple shooter wrapper

Author: Jonathan Richmond
C: 9/8/22
U: 8/30/23
"""

import LinearAlgebra
import MBD: MultipleShooter

export setPrintProgress!, solve!, solveUpdateEquation

"""
    setPrintProgress!(multipleShooter, print)

Return multiple shooter object with updated print progress

# Arguments
- `multipleShooter::MultipleShooter`: Multiple shooter object
- `print::Bool`: Print corrections progress?
"""
function setPrintProgress!(multipleShooter::MultipleShooter, print::Bool)
    multipleShooter.printProgress = print
end

"""
    solve!(multipleShooter, initialGuess)

Return converged solution

# Arguments
- `multipleShooter::MultipleShooter`: Multiple shooter object
- `initialGuess::MultipleShooterProblem`: Unsolved multiple shooter problem
"""
function solve!(multipleShooter::MultipleShooter, initialGuess::MBD.MultipleShooterProblem)
    buildProblem!(initialGuess)
    solutionInProgress::MBD.MultipleShooterProblem = deepClone(initialGuess)
    ((getNumberFreeVariables(solutionInProgress) == 0) || (getNumberConstraints(solutionInProgress) == 0)) ? (return solutionInProgress) : multipleShooter.recentIterationCount = 0
    while !isConverged(multipleShooter.convergenceCheck, solutionInProgress) && (multipleShooter.recentIterationCount < multipleShooter.maxIterations)
        if multipleShooter.recentIterationCount > 0
            freeVariableStep::Vector{Float64} = solveUpdateEquation(multipleShooter, solutionInProgress)
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
- `multipleShooter::MultipleShooter`: Multiple shooter object
- `multipleShooterProblem::MultipleShooterProblem`: Multiple shooter problem object
"""
function solveUpdateEquation(multipleShooter::MultipleShooter, multipleShooterProblem::MBD.MultipleShooterProblem)
    for generator::MBD.AbstractUpdateGenerator in multipleShooter.updateGenerators
        canGenerateUpdate(generator, multipleShooterProblem) && (return getFullUpdate(generator, multipleShooterProblem))
    end
    throw(ErrorException("No update generators were successful"))
end
