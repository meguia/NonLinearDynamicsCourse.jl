module NonLinearDynamicsCourse

    using Plots
    using LinearAlgebra
    using DifferentialEquations
    using ForwardDiff
    using IntervalRootFinding
    using StaticArrays 

    export  plot_nullclines,
            solve_plot_nullclines,
            solve_plot_nullclines_flux,   
            solve_plot_animated,
            classification_linear,
            plot_manifolds,
            phase_portrait,
            solve_plot_forced,
            poincare_forced,
            recurrence_plot

    include("DNL_utils.jl")

end
