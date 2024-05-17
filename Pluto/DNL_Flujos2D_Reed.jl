### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 27ea9f78-6276-4afd-bbd3-19587a4fd30e
using Pkg; Pkg.activate(joinpath(dirname(pwd()) , "src"))

# ╔═╡ 7c2fc40e-979f-4289-9b24-650a7364beac
push!(LOAD_PATH, joinpath(dirname(pwd()) , "src"))

# ╔═╡ 9e372d70-f6af-11ee-10b3-2150e19c755b
using DifferentialEquations, Plots, PlutoUI, JLD2, GLM, DataFrames

# ╔═╡ 33881071-dbb1-49a0-abf9-8ead17364c18
using NonLinearDynamicsCourse

# ╔═╡ f667d681-a2c4-4b04-a641-32692df1c5c2
function vreed!(du,u,p,t)
    (μ,k,v0) = p
	du[1] = u[2]
	du[2] =	(u[2] - v0)*(μ-(u[2] - v0)^2)-k*u[1]
end    

# ╔═╡ 9de8957e-2cab-4f33-a4e8-02932a704022
begin
	saved_values = load("freed.jld2")
	mu_v = saved_values["mu"]
	v0_v = saved_values["v0"]
	knotes = saved_values["knotes"]
	ampl_v = saved_values["ampl"]
	freq_v = saved_values["freq"]
end;	

# ╔═╡ c96c5652-386b-4129-8dc8-e76ef4769c26
begin 
	x = sqrt.(knotes[2:end])
	(nk, nv, nm) = size(freq_v)
	a0v = zeros(nv,nm)
	a1v = zeros(nv,nm)
	for n = 1:nv
		for m = 1:nm
			df1= DataFrame(x=x,y=freq_v[2:end,n,m])
			model1 = lm(@formula(y ~ x ), df1)
			(a0v[n,m], a1v[n,m]) = coef(model1)
		end
	end	
end	

# ╔═╡ b0e216b3-ba23-429b-bce3-2ad715313138
size(a0v)

# ╔═╡ e5d21f12-a369-493f-b83f-7d88235ff2dc
begin
	plot(mu_v,a0v',label="")
	plot!(mu_v,a0v[1,:],label="v0=0",lw=2,c=:black,xlabel="μ",ylabel="a0")
end	

# ╔═╡ 77007e68-52d2-4836-8732-2c758f64f061
begin
	plot(mu_v,a1v',label="")
	plot!(mu_v,a1v[1,:],label="v0=0",lw=2,c=:black,xlabel="μ",ylabel="a1")
end	

# ╔═╡ bc07de28-ce71-4b9e-85ab-b14ad8336e3c
begin
	plot(v0_v,a0v,label="")
	plot!(v0_v,a0v[:,161],label="mu="*string(mu_v[161]),lw=2,c=:orange)
	plot!(v0_v,a0v[:,131],label="mu="*string(mu_v[131]),lw=2,c=:cyan)
	plot!(v0_v,a0v[:,101],label="mu="*string(mu_v[101]),lw=2,c=:green)
	plot!(v0_v,a0v[:,71],label="mu="*string(mu_v[71]),lw=2,c=:red)
	plot!(v0_v,a0v[:,41],label="mu="*string(mu_v[41]),lw=2,c=:blue,xlabel="v0",ylabel="a0")
end	

# ╔═╡ 7c0ff248-cea7-49bd-b0b8-ff57ae860379
begin
	plot(v0_v,a1v,label="")
	plot!(v0_v,a1v[:,161],label="mu="*string(mu_v[161]),lw=2,c=:orange)
	plot!(v0_v,a1v[:,131],label="mu="*string(mu_v[131]),lw=2,c=:cyan)
	plot!(v0_v,a1v[:,101],label="mu="*string(mu_v[101]),lw=2,c=:green)
	plot!(v0_v,a1v[:,71],label="mu="*string(mu_v[71]),lw=2,c=:red)
	plot!(v0_v,a1v[:,41],label="mu="*string(mu_v[41]),lw=2,c=:blue,xlabel="v0",ylabel="a1")
end	

# ╔═╡ 13a1a5ca-5b80-4b89-b0f1-38f6e9ee44e6
f0 = 1.0

# ╔═╡ 68ee1490-50a7-40ea-87c4-bc7663fc5a79
kcomp = @. ((f0-a0v)/a1v)^2

# ╔═╡ 9da64564-73c9-4e42-b189-003a6c296afd
begin
	plot(v0_v,kcomp,label="",c=:gray)
	plot!(v0_v,kcomp[:,161],label="mu="*string(mu_v[161]),lw=2,c=:orange)
	plot!(v0_v,kcomp[:,131],label="mu="*string(mu_v[131]),lw=2,c=:cyan)
	plot!(v0_v,kcomp[:,101],label="mu="*string(mu_v[101]),lw=2,c=:green)
	plot!(v0_v,kcomp[:,71],label="mu="*string(mu_v[71]),lw=2,c=:red)
	plot!(v0_v,kcomp[:,41],label="mu="*string(mu_v[41]),lw=2,c=:blue)
	plot!(v0_v,kcomp[:,11],label="mu="*string(mu_v[11]),lw=2,c=:black,xlabel="v0",ylabel="k")
end	

# ╔═╡ dc0753f4-1add-4052-87a6-55c765d1a1f4


# ╔═╡ c0af1a99-7acf-4f05-8c98-52c3f898391d
md"""
x0 $(@bind x0_vreed Slider(-1.0:0.1:2.0,default=0.1;show_value=true)) \
y0 $(@bind y0_vreed Slider(-1.0:0.1:1.0,default=0.1;show_value=true)) \
μ : $(@bind μ_vreed Slider(-0.1:0.01:15.0,default=-0.1;show_value=true)) \
v0 : $(@bind v0_vreed Slider(-1.5:0.01:2.5,default=-0.1;show_value=true)) \
k : $(@bind k_vreed Slider(0.1:0.01:3.0,default=-0.1;show_value=true)) \
tmax : $(@bind tmax_vreed Slider(10:10:300.0,default=10.0;show_value=true)) 
"""

# ╔═╡ 85243d23-6d6f-4175-ab1b-a67d54309c10
begin
	sol3 = solve(ODEProblem(vreed!,[x0_vreed;y0_vreed],tmax_vreed,[μ_vreed , k_vreed, v0_vreed]),saveat=0.1)
	#plot(sol1, idxs=(0,1))
	plot(sol3, idxs=(0,2))
	
end	

# ╔═╡ 3c7f6407-bfc8-480f-922d-b98e5c5a4b69
flux2d_nullclines(vreed!,[x0_vreed;y0_vreed],tmax_vreed,[μ_vreed , k_vreed, v0_vreed],npts=41,xlims=[-6,6],ylims=[-1,2],title="Lengüeta (Rayleigh) Modificada")

# ╔═╡ 1dc81c7c-514a-4af0-a3c2-452ee9f71fd3
# ╠═╡ disabled = true
#=╠═╡
begin
	#esto es para generar la figura del paper
	sol1 = solve(ODEProblem(vreed!,[-0.1;0],tmax_vreed,[0.1 ,k_vreed, 0]),saveat=0.1)
	sol2 = solve(ODEProblem(vreed!,[-0.1;0],tmax_vreed,[0.79 ,k_vreed, 0]),saveat=0.1)
	plt1a = plot(sol1, idxs = (0,2),label="",ylabel="y",lc=:black,grid=false, title="μ=0.1, v₀=0")
	plt1b = plot(sol1, idxs = (1,2),label="",xlabel="x",ylabel="y",lc=:black,grid=false)
	plt2a = plot(sol2, idxs = (0,2),label="",ylabel="y",lc=:black,grid=false, title="μ=0.79, v₀=0")
	plt2b = plot(sol2, idxs = (1,2),label="",xlabel="x",ylabel="y",lc=:black,grid=false)
	plt3a = plot(sol3, idxs = (0,2),label="",ylabel="y",lc=:black,grid=false, title="μ=0.79, v₀=0.5")
	plt3b = plot(sol3, idxs = (1,2),label="",xlabel="x",ylabel="y",lc=:black,grid=false)
	plt_all = plot(plt1a,plt2a,plt3a,plt1b,plt2b,plt3b,layout=(2,3),size=(1200,700))
	savefig(plt_all, "fig2.png")
	plt_all
end	
  ╠═╡ =#

# ╔═╡ a0e34f60-e242-43af-87ae-17c9575ad1dc
html"""
<style>
main {
    max-width: 1200px;
}
input[type*="range"] {
	width: 90%;
}
</style>
"""

# ╔═╡ Cell order:
# ╠═9e372d70-f6af-11ee-10b3-2150e19c755b
# ╠═27ea9f78-6276-4afd-bbd3-19587a4fd30e
# ╠═33881071-dbb1-49a0-abf9-8ead17364c18
# ╠═7c2fc40e-979f-4289-9b24-650a7364beac
# ╠═f667d681-a2c4-4b04-a641-32692df1c5c2
# ╠═9de8957e-2cab-4f33-a4e8-02932a704022
# ╠═c96c5652-386b-4129-8dc8-e76ef4769c26
# ╠═b0e216b3-ba23-429b-bce3-2ad715313138
# ╠═e5d21f12-a369-493f-b83f-7d88235ff2dc
# ╠═77007e68-52d2-4836-8732-2c758f64f061
# ╠═bc07de28-ce71-4b9e-85ab-b14ad8336e3c
# ╠═7c0ff248-cea7-49bd-b0b8-ff57ae860379
# ╠═13a1a5ca-5b80-4b89-b0f1-38f6e9ee44e6
# ╠═68ee1490-50a7-40ea-87c4-bc7663fc5a79
# ╠═9da64564-73c9-4e42-b189-003a6c296afd
# ╠═dc0753f4-1add-4052-87a6-55c765d1a1f4
# ╠═85243d23-6d6f-4175-ab1b-a67d54309c10
# ╟─c0af1a99-7acf-4f05-8c98-52c3f898391d
# ╠═3c7f6407-bfc8-480f-922d-b98e5c5a4b69
# ╠═1dc81c7c-514a-4af0-a3c2-452ee9f71fd3
# ╠═a0e34f60-e242-43af-87ae-17c9575ad1dc
