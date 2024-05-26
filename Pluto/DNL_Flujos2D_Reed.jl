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

# ╔═╡ 146e013e-9448-4c89-9900-f5496fef8268
mu_v[80]

# ╔═╡ 2bb640a1-f234-4ce7-a14a-3e89ae5c8475
begin
	plot(v0_v,a0v[:,40]/2.3,label="a0",c=:red)
	plot!(v0_v,1 .-a1v[:,40],label="a1",c=:black,xlabel="v0")
end	

# ╔═╡ 4a414373-6810-46d1-a83a-51fc46a388d7
begin
	muidx = 16:5:60
	plot(v0_v,a0v[:,muidx],c=:red)
	plot!(v0_v,2.3 .*(1 .-a1v[:,muidx]),c=:black)
	scatter!(sqrt.(mu_v[muidx]/3),0 .*muidx)
	plot!(xlabel="v0",ylabel="a0",legend=false,xlims=(0,0.4))
end	

# ╔═╡ 204712ad-c798-4544-b16e-81619ed55f64
size(freq_v)

# ╔═╡ 13aeb37d-0f1d-4d83-b6dc-f720dc48d9e0
begin
	v0idx = 30
	kv = 0.25:0.01:4
	title0 = "Frecuencia en funcion de k para diferentes valores de mu con v0 = "*string(v0_v[v0idx])
	plt1 = plot(xlabel="k",ylabel="f",title=title0,size=(1000,400))
	for m in muidx
		plot!(plt1,kv, (a0v[v0idx,m].+a1v[v0idx,m]*sqrt.(kv)),c=:red,label="")
		plot!(plt1,kv, (a0v[v0idx,m].+(1.0-a0v[v0idx,m]/2.3)*sqrt.(kv)),c=:black,label="")
		scatter!(plt1,knotes,freq_v[:,v0idx,m],label=string(mu_v[m]),yaxis=log2,yticks=[0.5,1,2],xticks=[0.25,0.5,1,2,4],ylims=(0.45,2.2))
	end	
	plt1
end	

# ╔═╡ a2e4e348-5153-4dc2-8a09-34f87c0199c3
begin
	title1 = "Freq/sqrt(k) en funcion de k para diferentes valores de mu con v0 = "*string(v0_v[v0idx])
	plt2 = plot(xlabel="k",ylabel="f",title=title1,size=(1000,400))
	for m in muidx
		plot!(plt2,kv, (a0v[v0idx,m].+a1v[v0idx,m]*sqrt.(kv))./sqrt.(kv),c=:red,label="")
		plot!(plt2,kv, (a0v[v0idx,m].+(1.0-a0v[v0idx,m]/2.3)*sqrt.(kv))./sqrt.(kv),c=:black,label="")
		scatter!(knotes,freq_v[:,v0idx,m]./sqrt.(knotes),label=string(mu_v[m]))
	end	
	plt2
end	

# ╔═╡ ef00f387-2ef4-4e45-8687-cd5860e17f31
begin
	m2 = 31
	title2 ="Freq/sqrt(k) en funcion de k para diferentes valores de v0 con mu = "*string(mu_v[m2])
	plt3 = plot(xlabel="k",ylabel="f",title=title2,size=(1000,400))
	for v2 in 1:10:51
		plot!(plt3,kv, (a0v[v2,m2].+a1v[v2,m2]*sqrt.(kv))./sqrt.(kv),c=:red,label="")
		plot!(plt3,kv, (a0v[v2,m2].+(1.0-a0v[v2,m2]/2.3)*sqrt.(kv))./sqrt.(kv),c=:black,label="")
		scatter!(plt3,knotes,freq_v[:,v2,m2]./sqrt.(knotes),label=string(v0_v[v2]))
	end	
	plt3
end	

# ╔═╡ c0af1a99-7acf-4f05-8c98-52c3f898391d
md"""
x0 $(@bind x0_vreed Slider(-1.0:0.1:2.0,default=0.1;show_value=true)) \
y0 $(@bind y0_vreed Slider(-1.0:0.1:1.0,default=0.1;show_value=true)) \
μ : $(@bind μ_vreed Slider(-0.1:0.01:15.0,default=-0.1;show_value=true)) \
v0 : $(@bind v0_vreed Slider(-1.5:0.01:2.5,default=-0.1;show_value=true)) \
k : $(@bind k_vreed Slider(0.1:0.01:3.0,default=-0.1;show_value=true)) \
tmax : $(@bind tmax_vreed Slider(10:10:300.0,default=10.0;show_value=true)) 
"""

# ╔═╡ 1dc81c7c-514a-4af0-a3c2-452ee9f71fd3
begin
	#esto es para generar la figura del paper
	k = 0.28
	μ = 0.8
	v0 = 0.45
	tmax = 100
	sol1 = solve(ODEProblem(vreed!,[-0.1;0],tmax,[0.1 ,k, 0]),saveat=0.1)
	sol2 = solve(ODEProblem(vreed!,[-0.1;0],tmax,[μ ,k, 0]),saveat=0.1)
	sol3 = solve(ODEProblem(vreed!,[-0.1;0],tmax,[μ , k, v0]),saveat=0.1)
	plt1a = plot(sol1, idxs = (0,2),label="",ylabel="y",lc=:black,grid=false, title="μ=0.1, v₀=0")
	plt1b = plot(sol1, idxs = (1,2),label="",xlabel="x",ylabel="y",lc=:black,grid=false)
	plt2a = plot(sol2, idxs = (0,2),label="",ylabel="y",lc=:black,grid=false, title="μ=$μ, v₀=0")
	plt2b = plot(sol2, idxs = (1,2),label="",xlabel="x",ylabel="y",lc=:black,grid=false)
	plt3a = plot(sol3, idxs = (0,2),label="",ylabel="y",lc=:black,grid=false, title="μ=$μ, v₀=$v0")
	plt3b = plot(sol3, idxs = (1,2),label="",xlabel="x",ylabel="y",lc=:black,grid=false)
	plt_all = plot(plt1a,plt2a,plt3a,plt1b,plt2b,plt3b,layout=(2,3),size=(1200,700))
	savefig(plt_all, "fig2.png")
	plt_all
end	

# ╔═╡ 88e254d0-60c3-4438-b6e8-37a1c334cead
knotes

# ╔═╡ 6a233c53-92a5-42b8-a113-77893cb42df0
begin
	ik = findmin(abs.(knotes .- k))[2]
	fig3 = contourf(v0_v,mu_v,freq_v[ik,:,:]',levels=levels=freq_v[ik,1,1]*2 .^(-12/12:1/12:0),colorbar=false)
	plot!(fig3,v0_v,3*v0_v.^2,c=:black,lw=3,label="Hopf")
	scatter!(fig3,[v0],[μ],c=:blue,label="")
	scatter!(fig3,[0],[μ],c=:blue,label="")
	scatter!(fig3,[0],[0.1],c=:blue,label="",xlabel="v₀",ylabel="μ")
	savefig(fig3, "fig3.png")
	fig3
end	

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

# ╔═╡ 0064e75e-ed6b-4fed-a854-326d20edd543
size(freq_v)

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
# ╠═146e013e-9448-4c89-9900-f5496fef8268
# ╠═2bb640a1-f234-4ce7-a14a-3e89ae5c8475
# ╠═4a414373-6810-46d1-a83a-51fc46a388d7
# ╠═204712ad-c798-4544-b16e-81619ed55f64
# ╠═13aeb37d-0f1d-4d83-b6dc-f720dc48d9e0
# ╠═a2e4e348-5153-4dc2-8a09-34f87c0199c3
# ╠═ef00f387-2ef4-4e45-8687-cd5860e17f31
# ╟─c0af1a99-7acf-4f05-8c98-52c3f898391d
# ╠═1dc81c7c-514a-4af0-a3c2-452ee9f71fd3
# ╠═88e254d0-60c3-4438-b6e8-37a1c334cead
# ╠═6a233c53-92a5-42b8-a113-77893cb42df0
# ╟─a0e34f60-e242-43af-87ae-17c9575ad1dc
# ╠═0064e75e-ed6b-4fed-a854-326d20edd543
