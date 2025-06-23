"""
Homotopy parameter end check wrapper

Author: Jonathan Richmond
C: 6/23/25
"""

import Logging
import MBD: HomotopyEndCheck

export isContinuationDone

"""
    isContinuationDone(homotopyCheck, data)

Return true if continuation is done

# Arguments
- `homotopyCheck::HomotopyEndCheck`: Homotopy parameter end check object
- `data::BCR4BP12ContinuationData`: BCR4BP P1-P2 continuation data object
"""
function isContinuationDone(homotopyCheck::HomotopyEndCheck, data::MBD.BCR4BP12ContinuationData)
    systemData::MBD.BCR4BPSystemData = data.previousSolution.nodes[1].dynamicsModel.systemData
    currentParam::Float64 = systemData.P4Mass/systemData.primaryData[3].mass
    maxParam::Float64 = homotopyCheck.maxParam

    Logging.@debug "Checking if homotopy parameter is reached"

    if currentParam >= maxParam
        Logging.@info "Homotopy parameter reached: $currentParam ≥ $maxParam"
        println("Homotopy parameter reached!")

        return true
    end

    Logging.@debug "Continuation ongoing: homotopy parameter $currentParam < $maxParam"
    return false
end
