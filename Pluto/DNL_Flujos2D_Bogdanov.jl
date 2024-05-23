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

# ╔═╡ 89db1067-13ac-48e2-a260-c1b4fe154c0f
Pkg.add("PolynomialRoots")

# ╔═╡ 7c2fc40e-979f-4289-9b24-650a7364beac
push!(LOAD_PATH, joinpath(dirname(pwd()) , "src"))

# ╔═╡ 9e372d70-f6af-11ee-10b3-2150e19c755b
using DifferentialEquations, Plots, PlutoUI, JLD2, GLM, DataFrames, PolynomialRoots

# ╔═╡ 33881071-dbb1-49a0-abf9-8ead17364c18
using NonLinearDynamicsCourse

# ╔═╡ 6bc1183a-9f68-4f1d-bcf8-37e7ae15fc67
using Random

# ╔═╡ f667d681-a2c4-4b04-a641-32692df1c5c2
function bogdanov!(du, u, p, t)
	(μ1,μ2,δ, η) = p 
	du[1] = u[2]
    du[2] = μ1+u[1]*(μ2-δ*u[2]+u[1]*(1-u[1]-u[2]))
end

# ╔═╡ 7ea8fe01-24ff-4fd3-83d7-b05ea831b268
function noise!(du,u,p,t)
	du[1] = p[end]
	du[2] = 0
end

# ╔═╡ 9de8957e-2cab-4f33-a4e8-02932a704022
begin
	file_temp = load("Homall.jld2")
	dval = file_temp["dvalues"]
	hom = file_temp["hom"]
end;	

# ╔═╡ 524eebfa-1bea-446a-a877-de1f0b42b313
dval

# ╔═╡ 12d3ffee-818e-4411-8769-eef14bc3be4b
begin
	sn1(μ2) =	-(3*μ2+2/3+sqrt(3*μ2+1.0)*(2*μ2+2/3))/9.0
	sn2(μ2) =	-(3*μ2+2/3-sqrt(3*μ2+1.0)*(2*μ2+2/3))/9.0
	h2(μ2,δ) = δ*(μ2-δ-δ^2)
	function bt(δ) 
		δ2 = 9*(δ+3*δ^2+3*δ^3)
		m2 = (real(roots([-(δ2+1),0,3.0+9*δ,-2.0])).^2 .-1 )./3
		m2[1]
	end	
end;	

# ╔═╡ 77ec9c72-454f-4c82-9cb9-fdc97803f31a
md"""
x0 $(@bind x0 Slider(-1:0.1:1,default=0.1;show_value=true)) \
y0 $(@bind y0 Slider(-1:0.1:1,default=0.1;show_value=true)) \
dμ2 : $(@bind dμ2 Slider(-0.02:0.001:0.02,default=0.02;show_value=true))\
μ20 : $(@bind μ20 Slider(-0.33:0.001:0.4,default=-0.02;show_value=true)) \
δ : $(@bind δ Slider(0.0:0.01:2.0,default=0.5;show_value=true)) \
η : $(@bind η Slider(0.0:0.001:0.1,default=0.0;show_value=true)) \
tmax $(@bind tmax Slider(30:10:1000,default=100;show_value=true)) \
"""

# ╔═╡ e8bd6226-52de-47c3-af77-b672ff7e256f
begin
	bt2 = bt(δ)
	μsn = -0.33:0.01:max(bt2,0.8)
	μh1 = -0.5:0.01:0
	μh2 = -0.5:0.01:0.8
	nhom = argmin(abs.(δ .- dval/10))
	plta = plot(μsn,sn1.(μsn),c=:black,label="Saddle Node",xlabel="μ2",ylabel="μ1",title="δ=1.5")
	plot!(plta,μsn,sn2.(μsn),c=:black,label="")
	plot!(plta,μh1,μh1*0,c=:black,ls=:dash,label="Hopf")
	#plot!(plta,μh2,h2.(μh2,δ),c=:black,ls=:dash,label="")
	plot!(plta,hom[12][2,:],hom[12][1,:],c=:black,ls=:dot,label="Homoclinic")
	scatter!(plta,[0],[0],c=:gray,ms=3,label="BT")
	scatter!(plta,[bt2],[sn2(bt2)],c=:gray,ms=3,label="")
	#scatter!([μ20-dμ2],[sn1(μ20)],c=:blue,label="",xlims=(-0.4,0.3),ylims=(-0.3,0.05))
	scatter!(plta,[0,-0.3],[-0.15,-0.005],ms=3,c=:blue,label="η=0")
	scatter!(plta,[0,-0.3],[-0.145,0.005],ms=3,c=:red,label="η=0.02",xlims=(-0.4,0.3),ylims=(-0.3,0.05))
	pltb = plot(μsn,sn1.(μsn),c=:black,label="Saddle Node",xlabel="μ2",ylabel="μ1",title="δ=0.5")
	plot!(pltb,μsn,sn2.(μsn),c=:black,label="")
	plot!(pltb,μh1,μh1*0,c=:black,ls=:dash,label="Hopf")
	plot!(pltb,μh2,h2.(μh2,δ),c=:black,ls=:dash,label="")
	plot!(pltb,hom[5][2,:],hom[5][1,:],c=:black,ls=:dot,label="Homoclinic")
	scatter!(pltb,[0],[0],c=:gray,ms=3,label="BT")
	scatter!(pltb,[bt2],[sn2(bt2)],c=:gray,ms=3,label="")
	#scatter!([μ20-dμ2],[sn1(μ20)],c=:blue,label="",xlims=(-0.4,0.3),ylims=(-0.3,0.05))
	scatter!(pltb,[0],[-0.15],ms=3,c=:blue,label="η=0")
	scatter!(pltb,[0],[-0.145],ms=3,c=:red,label="η=0.02",xlims=(-0.4,0.3),ylims=(-0.3,0.05))
	plt_ab = plot(pltb,plta,layout=(1,2),size=(1200,500))
	savefig(plt_ab, "fig6.png")
	plt_ab
end	

# ╔═╡ 85243d23-6d6f-4175-ab1b-a67d54309c10
begin
	# Bogdanov without noise
	sol1 = solve(SDEProblem(bogdanov!,noise!,[-0.2;0],220,[-0.15,0,1.5,0]),saveat=0.1)
	plt2a = plot(sol1, idxs=(0,2),label="",ylabel="y",lc=:black,grid=false, title="μ1=-0.15, μ2=0, δ=1.5")
	plt2b = plot(sol1, idxs = (1,2),label="",xlabel="x",ylabel="y",lc=:black,grid=false)
	sol2 = solve(SDEProblem(bogdanov!,noise!,[-0.3;0],220,[-0.15,0,0.5,0]),saveat=0.1)
	plt1a = plot(sol2, idxs=(0,2),label="",ylabel="y",lc=:black,grid=false, title="μ1=-0.15, μ2=0, δ=0.5")
	plt1b = plot(sol2, idxs = (1,2),label="",xlabel="x",ylabel="y",lc=:black,grid=false)
	sol3 = solve(SDEProblem(bogdanov!,noise!,[-0.2;0],220,[-0.005,-0.3,1.5,0]),saveat=0.1)
	plt3a = plot(sol3, idxs=(0,2),label="",ylabel="y",lc=:black,grid=false, title="μ1=-0.005, μ2=-0.3, δ=1.5")
	plt3b = plot(sol3, idxs = (1,2),label="",xlabel="x",ylabel="y",lc=:black,grid=false)
	plt_all = plot(plt1a,plt2a,plt3a,plt1b,plt2b,plt3b,layout=(2,3),size=(1200,700))
	savefig(plt_all, "fig4.png")
	plt_all
end	

# ╔═╡ 7c660bce-d604-4d16-8c19-614fd3db46e9
begin
	# Bogdanov with noise
	Random.seed!(1234)
	sol4 = solve(SDEProblem(bogdanov!,noise!,[0.7;0],500,[-0.145,0.0,1.5,0.02]),saveat=0.02)
	plt5a = plot(sol4, idxs=(0,2),label="",ylabel="y",lc=:black,grid=false, title="μ1=-0.145, μ2=0, δ=1.5")
	plt5b = plot(sol4, idxs = (1,2),label="",xlabel="x",ylabel="y",lc=:black,grid=false)
	sol5 = solve(SDEProblem(bogdanov!,noise!,[0.65;0.0],500,[-0.145,0.0,0.5,0.02]),saveat=0.02)
	plt4a = plot(sol5, idxs=(0,2),label="",ylabel="y",lc=:black,grid=false, title="μ1=-0.145, μ2=0, δ=0.5")
	plt4b = plot(sol5, idxs = (1,2),label="",xlabel="x",ylabel="y",lc=:black,grid=false)
	sol6 = solve(SDEProblem(bogdanov!,noise!,[0.2;0],500,[0.005,-0.3,1.5,0.02]),saveat=0.1)
	plt6a = plot(sol6, idxs=(0,2),label="",ylabel="y",lc=:black,grid=false, title="μ1=0.005, μ2=-0.3, δ=1.5")
	plt6b = plot(sol6, idxs = (1,2),label="",xlabel="x",ylabel="y",lc=:black,grid=false)
	plt_all2 = plot(plt4a,plt5a,plt6a,plt4b,plt5b,plt6b,layout=(2,3),size=(1200,700))
	savefig(plt_all2, "fig5.png")
	plt_all2
end	

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
# ╠═89db1067-13ac-48e2-a260-c1b4fe154c0f
# ╠═9e372d70-f6af-11ee-10b3-2150e19c755b
# ╠═27ea9f78-6276-4afd-bbd3-19587a4fd30e
# ╠═33881071-dbb1-49a0-abf9-8ead17364c18
# ╠═7c2fc40e-979f-4289-9b24-650a7364beac
# ╠═f667d681-a2c4-4b04-a641-32692df1c5c2
# ╠═7ea8fe01-24ff-4fd3-83d7-b05ea831b268
# ╠═9de8957e-2cab-4f33-a4e8-02932a704022
# ╠═524eebfa-1bea-446a-a877-de1f0b42b313
# ╠═12d3ffee-818e-4411-8769-eef14bc3be4b
# ╠═e8bd6226-52de-47c3-af77-b672ff7e256f
# ╟─77ec9c72-454f-4c82-9cb9-fdc97803f31a
# ╠═85243d23-6d6f-4175-ab1b-a67d54309c10
# ╠═6bc1183a-9f68-4f1d-bcf8-37e7ae15fc67
# ╠═7c660bce-d604-4d16-8c19-614fd3db46e9
# ╠═1dc81c7c-514a-4af0-a3c2-452ee9f71fd3
# ╠═a0e34f60-e242-43af-87ae-17c9575ad1dc
