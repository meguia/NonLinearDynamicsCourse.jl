### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 2050b25a-537c-4057-a930-4d8f713ab9d2
using Pkg; Pkg.add("HclinicBifurcationKit")

# ╔═╡ 9e372d70-f6af-11ee-10b3-2150e19c755b
using DifferentialEquations, Plots, Parameters, Setfield, BifurcationKit, PlutoUI

# ╔═╡ 2a9addce-3887-4ae7-903f-697bc1ff0dec
function lff(u, p)
	@unpack ϵ1,ϵ2 = p 
	[
		u[2]
		ϵ1-u[2]+u[1]*(1.0+ϵ2*u[1]+u[2]-u[1]^2)
	]	
end

# ╔═╡ 71656988-ec80-460d-b887-fc2551e34529
par = (ϵ1 = 0.1, ϵ2 = 0.3)

# ╔═╡ 3fdb9050-22b9-4d42-8785-1fb96d08f117
prob = BifurcationProblem(lff, [0.1; 0.0], par, (@lens _.ϵ1),record_from_solution = (x, p) -> x[1])

# ╔═╡ e48feb36-124e-41bf-a65f-871b852ca23a
opts_br = ContinuationPar(p_max = 0.5, p_min = 0.1,dsmax=0.01,dsmin=1e-4,ds=1e-4,max_steps = 100,n_inversion = 8);

# ╔═╡ e8c305eb-01c3-4f5f-a227-0c716edba03f
br = continuation(prob, PALC(), opts_br, bothside = true)

# ╔═╡ fa0eac04-8709-484f-87e1-ae18cea1ccba
opts_br2 = ContinuationPar(p_max = 1.0, p_min = 0.0,dsmax=0.01,dsmin=1e-4,ds=1e-4,max_steps = 100,n_inversion = 4, detect_event = 0,detect_bifurcation=3);

# ╔═╡ c1296cb7-9acd-44d4-bff2-9dc411ef0317
snlc = continuation(br, 2, (@lens _.ϵ2), opts_br2, normC = norminf, bothside = true, update_minaug_every_step = 2, detect_codim2_bifurcation = 2)

# ╔═╡ d4fd2dd5-c3bd-49d5-887b-497821c04809
plot(snlc)

# ╔═╡ 4eea2bee-8de4-4679-bc72-aeb65c56203d
md"""
# Takens Bogdanov
"""

# ╔═╡ a0e34f60-e242-43af-87ae-17c9575ad1dc
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
# ╠═9e372d70-f6af-11ee-10b3-2150e19c755b
# ╠═2050b25a-537c-4057-a930-4d8f713ab9d2
# ╠═2a9addce-3887-4ae7-903f-697bc1ff0dec
# ╠═71656988-ec80-460d-b887-fc2551e34529
# ╠═3fdb9050-22b9-4d42-8785-1fb96d08f117
# ╠═e48feb36-124e-41bf-a65f-871b852ca23a
# ╠═e8c305eb-01c3-4f5f-a227-0c716edba03f
# ╠═fa0eac04-8709-484f-87e1-ae18cea1ccba
# ╠═c1296cb7-9acd-44d4-bff2-9dc411ef0317
# ╠═d4fd2dd5-c3bd-49d5-887b-497821c04809
# ╟─4eea2bee-8de4-4679-bc72-aeb65c56203d
# ╟─a0e34f60-e242-43af-87ae-17c9575ad1dc
