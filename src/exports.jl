"""
Package exports

Author: Jonathan LeFevre Richmond
C: 4/14/26
U: 8/7/26
"""


export
    # Types
    BodyData,
    SystemData,
    AbstractDynamicsModel, CR3BPDynamicsModel,
    # AbstractEquationsOfMotion, CR3BPEquationsOfMotion,

    # Dynamics
    # SystemData methods
    getNumPrimaries,
    # DynamicsModel methods
    adjustInitialConditions, appendExtraInitialConditions,
    extractStateTransitionMatrix, getCharLengths, getCharMasses, getCharTimes,
    getEquilibriumPoint, getExcursion, getHamiltonian, getJacobiConstant,
    getLinearVariation, getMassRatios, getParameterDependencies,
    getPrimaryState, getPseudopotential, getPseudopotentialHessian,
    getPseudopotentialJacobian, getStateSize, getTidalAcceleration,
    # getEpochDependencies, getEquationsOfMotion, isEpochIndependent

    # Utilities
    # SPICE methods
    getIDCode
