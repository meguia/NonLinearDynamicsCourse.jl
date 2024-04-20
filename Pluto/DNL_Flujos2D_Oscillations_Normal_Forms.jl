### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 878ca803-f232-45d4-8f2a-e72fd14c0324
using Pkg; Pkg.add("PyPlot")

# ╔═╡ 72fa8994-47dc-4e0c-80fd-021f5fc629e9
using DifferentialEquations, Plots

# ╔═╡ 163d352d-ce42-4828-a418-520033c23719
# Escribimos la ecuacion de la lengueta
function vreed!(du,u,p,t)
    (μ,k,v0) = p
	v = u[2] - v0
    du[1] = u[2]
    du[2] = v*(μ-v*v)-k*u[1]
    du
end    

# ╔═╡ fca48a5c-d1f0-4069-85c9-71c6b612b25f
# Para la Figura 1
begin
	k = 2.0
	μ = 0.4
	sol1a = solve(ODEProblem(vreed!,[0.0;1.0],(0.0,5.0),[-0.4,k,0.0]));
	sol1b = solve(ODEProblem(vreed!,[0.0;1.0],(0.0,20.0),[-0.4,k,0.0]));
	sol2a = solve(ODEProblem(vreed!,[0.02;0.0],(0.0,20.0),[μ,k,0.0]));
	sol2b = solve(ODEProblem(vreed!,[0.02;0.0],(0.0,14.0),[μ,k,0.0]));
	sol2c = solve(ODEProblem(vreed!,[0.0;1.0],(0.0,5.8),[μ,k,0.0]));
	sol2d = solve(ODEProblem(vreed!,[0.0;0.727],(0.0,8.5),[μ,k,0.0]));
end;

# ╔═╡ e58881c1-3f0c-4b82-b787-37325442226e
begin
	xlims = (-1,1)
	ylims = (-1,1)
	plt1 = scatter([0],[0],color=:black,xlims=xlims,ylims=ylims,legend=false)
	plot!(plt1,sol1a,idxs=(1,2),c=:black,arrow=true,ls=:dash,lw=1)
	plot!(plt1,sol1b,idxs=(1,2),c=:black,arrow=false,ls=:dash,lw=1)
	plot!(plt1,xlabel="",ylabel="",border=:none)
	plt2 = scatter([0],[0],color=:white,xlims=xlims,ylims=ylims,legend=false)
	plot!(plt2,sol2a,idxs=(1,2),c=:black,arrow=true,ls=:dash,lw=1)
	plot!(plt2,sol2b,idxs=(1,2),c=:black,arrow=true,ls=:dash,lw=1)
	plot!(plt2,sol2c,idxs=(1,2),c=:black,arrow=true,ls=:dash,lw=1)
	plot!(plt2,sol2d,idxs=(1,2),c=:black,arrow=true,lw=2)
	plot!(plt2,xlabel="",ylabel="",border=:none)
	plt = plot(plt1,plt2,layout = (1, 2),size=(1000,400))
	savefig(plt, "fig1.eps")
end	

# ╔═╡ 4f04d9a1-3e58-4ead-b5c4-daf4d88e9f29
html"""
<style>
main {
    max-width: 1000px;
}
input[type*="range"] {
	width: 90%;
}
</style>
"""

# ╔═╡ Cell order:
# ╠═878ca803-f232-45d4-8f2a-e72fd14c0324
# ╠═72fa8994-47dc-4e0c-80fd-021f5fc629e9
# ╠═163d352d-ce42-4828-a418-520033c23719
# ╠═fca48a5c-d1f0-4069-85c9-71c6b612b25f
# ╠═e58881c1-3f0c-4b82-b787-37325442226e
# ╟─4f04d9a1-3e58-4ead-b5c4-daf4d88e9f29
