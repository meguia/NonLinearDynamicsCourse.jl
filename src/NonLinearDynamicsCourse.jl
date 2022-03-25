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
            attractor_basin,
            solve_plot_forced,
            poincare_forced,
            poincare_forced_zoom,
            recurrence_plot,
            saddle_orbit2D,
            saddle_manifolds_forced

    include("DNL_utils.jl")

end
