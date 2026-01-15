"""
Core types

Author: Jonathan LeFevre Richmond
C: 12/12/25
U: 12/24/25
"""


"""Enumerated types"""
@enum EquationType begin
    ARCLENGTH
    FULL
    MOMENTUM
    SIMPLE
    STM
end

@enum ModelType begin
    CR3BP
end
