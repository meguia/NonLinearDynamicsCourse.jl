### A Pluto.jl notebook ###
# v0.19.4

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

# ╔═╡ 7c860884-5afd-4b4a-b358-4f2bce8bd395
using Pkg; Pkg.activate("../..")

# ╔═╡ 38c4e220-d708-11ec-3968-fbbd41c26155
using PlutoUI, Plots, DifferentialEquations, BifurcationKit, Setfield

# ╔═╡ ee3d45e6-da8d-433d-a360-ef2d2bc353b6
using NonLinearDynamicsCourse

# ╔═╡ 62cba86d-4406-4d16-8715-b29bee7d3178
md"""
# Logistic Equation

$\dot{x}=Rx\left(1-\frac{x}{K}\right)$
"""

# ╔═╡ 4dc32b5a-5794-4d2d-844e-cdb80bec2e7b
# Logistic Equation
logistic(x,p,t)=p[1]*x*(1.0-x/p[2])

# ╔═╡ 77766c99-4142-4acb-a855-ab133e67aa66
@bind pars2 (
	PlutoUI.combine() do bind
		md"""
		R: $(bind(Slider(0:0.02:2.0,default=0.1;show_value=true))) \
		K: $(bind(Slider(0.01:0.02:2.0,default=0.1;show_value=true))) \
		x0: $(bind(Slider(0:0.02:2.0,default=0.1;show_value=true)))
		"""
	end
)

# ╔═╡ 16a156b5-92ab-4e1c-ab70-c3b348186c57
flux1D(logistic,pars2[3],100.0,pars2;xlims=[-0.1,2.0],title="Logistic")

# ╔═╡ ce4da9bf-316b-45d2-8504-f2a63cdda247
md"""
# Logistic Equation with Harvest

$\dot{x}=Rx\left(1-\frac{x}{K}\right)-H$
"""

# ╔═╡ 29474ca0-839b-45de-a7bf-7d91c5ea59aa
logharvest1(x,p,t)=p[1]*x*(1.0-x/p[2])-p[3]

# ╔═╡ 8fc407f4-a242-4cdc-b7e0-062c8f5782d1
@bind pars_harvest (
	PlutoUI.combine() do bind
		md"""
		R: $(bind(Slider(0:0.02:2.0,default=0.1;show_value=true))) \
		K: $(bind(Slider(0.01:0.02:2.0,default=0.1;show_value=true))) \
		H: $(bind(Slider(0.01:0.02:2.0,default=0.1;show_value=true))) \
		x0: $(bind(Slider(0:0.02:2.0,default=0.1;show_value=true)))
		"""
	end
)

# ╔═╡ 794936ae-d9ae-4319-994f-4282c99c4238
flux1D(logharvest1,pars_harvest[4],300.0,pars_harvest,(u)->(u<0);xlims=[0.0,2.0],title="Logistic with Harvest")

# ╔═╡ Cell order:
# ╠═7c860884-5afd-4b4a-b358-4f2bce8bd395
# ╠═38c4e220-d708-11ec-3968-fbbd41c26155
# ╠═ee3d45e6-da8d-433d-a360-ef2d2bc353b6
# ╟─62cba86d-4406-4d16-8715-b29bee7d3178
# ╠═4dc32b5a-5794-4d2d-844e-cdb80bec2e7b
# ╟─16a156b5-92ab-4e1c-ab70-c3b348186c57
# ╟─77766c99-4142-4acb-a855-ab133e67aa66
# ╟─ce4da9bf-316b-45d2-8504-f2a63cdda247
# ╠═29474ca0-839b-45de-a7bf-7d91c5ea59aa
# ╟─794936ae-d9ae-4319-994f-4282c99c4238
# ╟─8fc407f4-a242-4cdc-b7e0-062c8f5782d1
