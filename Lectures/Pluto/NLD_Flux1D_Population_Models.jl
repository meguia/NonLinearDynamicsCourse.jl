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
using PlutoUI, Plots, DifferentialEquations, NonLinearDynamicsCourse

# ╔═╡ 6b6a5e19-5298-4969-b1c9-699aa1cb2996
html"""
<style>
input[type*="range"] {
	width: 50%;
}
</style>
"""

# ╔═╡ 62cba86d-4406-4d16-8715-b29bee7d3178
md"""
# Logistic Equation

$\dot{x}=Rx\left(1-\frac{x}{K}\right)$
"""

# ╔═╡ 6d64cc27-1059-4dd7-8f17-1c6fd9a2eca0
html"""
<div style="position: relative; right: 0; top: 0; z-index: 300;"><iframe src="https://www.youtube.com/embed/JwYhhnuuINk" width=500 height=250  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></div>
"""

# ╔═╡ 4dc32b5a-5794-4d2d-844e-cdb80bec2e7b
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

# ╔═╡ b13d261e-1c1b-41b0-ae86-e20cd2304727
md"""
## Critical Slowing Down
"""

# ╔═╡ c3757398-060c-412c-997b-58331b2dee04
@bind pars_csd (
	PlutoUI.combine() do bind
		md"""
		H: $(bind(Slider(0.2:0.002:0.25,default=0.1;show_value=true))) \
		S: $(bind(Slider(0:0.001:0.1,default=0.0;show_value=true)))
		"""
	end
)

# ╔═╡ f1f04dd4-0447-4e48-9bab-f34058a6a0f1
flux1D(logharvest1,0.5+sqrt(0.25-pars_csd[1]),200.0,[1.0,1.0,pars_csd[1]],10.0,pars_csd[2],(u)->(u<0);xlims=[0.0,1.0],title="Log whith Harvest perturbed")

# ╔═╡ c803ed90-958c-4d0e-84d5-3e8cc9c8fbfe
md"""
# Consumer Equation

$\dot{x} = Rx\left(1-\frac{x}{K}\right) - Px$ 
"""

# ╔═╡ b4ea8364-a582-443c-bdb1-46e75a5c4d4e
# Consumer Equation
consumer(x,p,t)=p[1]*x*(1.0-x/p[2])-p[3]*x

# ╔═╡ 5afa3aff-af15-429b-a0fb-9a0dfcad740d
@bind pars_consumer (
	PlutoUI.combine() do bind
		md"""
		R: $(bind(Slider(0:0.02:1.0,default=0.1;show_value=true))) \
		K: $(bind(Slider(0.01:0.02:2.0,default=0.1;show_value=true))) \
		P: $(bind(Slider(0.0:0.01:0.5,default=0.1;show_value=true))) \
		x0: $(bind(Slider(0:0.02:2.0,default=0.1;show_value=true)))
		"""
	end
)

# ╔═╡ 5018d342-0a51-4275-8c0f-1b5d253c2ec0
flux1D(consumer,pars_consumer[4],300.0,pars_consumer,xlims=[0.0,2.0],title="Consumer Equation")

# ╔═╡ 84a596c1-3a0a-4edd-a7ea-f3fcd8839c2a
md"""
# Logistic Equation with Outbreak

$\dot{x} = Rx\left(1-\displaystyle\frac{x}{K}\right)-P\displaystyle\frac{x^2}{1+x^2}$
"""

# ╔═╡ 91b2718a-a921-4bae-99c0-91ec9c7e6479
html"""
<div style="position: relative; right: 0; top: 0; z-index: 300;"><iframe src="https://www.youtube.com/embed/1CSKTCS6st8" width=500 height=250  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></div>
"""

# ╔═╡ 5a0b156e-b0b2-43b8-b15e-bcdb65c154bf
logoutbreak(x,p,t)=p[1]*x*(1.0-x/p[2])-p[3]*x*x/(1+x*x)

# ╔═╡ 6528501f-3bb2-48a7-b24e-d32400033538
@bind pars_outbreak (
	PlutoUI.combine() do bind
		md"""
		R: $(bind(Slider(0:0.02:2.0,default=0.1;show_value=true))) \
		K: $(bind(Slider(0.01:0.02:10.0,default=0.1;show_value=true))) \
		P: $(bind(Slider(0.0:0.01:1.0,default=0.1;show_value=true))) \
		x0: $(bind(Slider(0:0.02:8.0,default=0.1;show_value=true)))
		"""
	end
)

# ╔═╡ a3270d8b-f167-4fb9-bcc5-c62c91b39649
flux1D(logoutbreak,pars_outbreak[4],300.0,pars_outbreak;xlims=[-0.2,8.0],title="Log with Outbreak")

# ╔═╡ 0d0d6d2f-56b2-4e01-b9a4-6e469eac90ad
html"""
<div style="position: relative; right: 0; top: 0; z-index: 300;"><iframe src="https://www.youtube.com/embed/MlwAI3BlDsU" width=500 height=250  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></div>
"""

# ╔═╡ 76cfd5c8-6298-4ae2-9d2d-9e60deb12f49
html"""
<div style="position: relative; right: 0; top: 0; z-index: 300;"><iframe src="https://www.youtube.com/embed/yUEXpUyi404" width=500 height=250  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></div>
"""

# ╔═╡ 0bf40552-b627-4fa7-b3ba-1d8c7ed4a568
html"""
<div style="position: relative; right: 0; top: 0; z-index: 250;"><iframe src="https://www.youtube.com/embed/Vu9oNWXv4Uk" width=500 height=250  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></div>
"""

# ╔═╡ Cell order:
# ╟─7c860884-5afd-4b4a-b358-4f2bce8bd395
# ╠═38c4e220-d708-11ec-3968-fbbd41c26155
# ╟─6b6a5e19-5298-4969-b1c9-699aa1cb2996
# ╟─62cba86d-4406-4d16-8715-b29bee7d3178
# ╟─6d64cc27-1059-4dd7-8f17-1c6fd9a2eca0
# ╠═4dc32b5a-5794-4d2d-844e-cdb80bec2e7b
# ╟─16a156b5-92ab-4e1c-ab70-c3b348186c57
# ╟─77766c99-4142-4acb-a855-ab133e67aa66
# ╟─ce4da9bf-316b-45d2-8504-f2a63cdda247
# ╠═29474ca0-839b-45de-a7bf-7d91c5ea59aa
# ╟─794936ae-d9ae-4319-994f-4282c99c4238
# ╟─8fc407f4-a242-4cdc-b7e0-062c8f5782d1
# ╟─b13d261e-1c1b-41b0-ae86-e20cd2304727
# ╟─f1f04dd4-0447-4e48-9bab-f34058a6a0f1
# ╟─c3757398-060c-412c-997b-58331b2dee04
# ╟─c803ed90-958c-4d0e-84d5-3e8cc9c8fbfe
# ╠═b4ea8364-a582-443c-bdb1-46e75a5c4d4e
# ╟─5018d342-0a51-4275-8c0f-1b5d253c2ec0
# ╟─5afa3aff-af15-429b-a0fb-9a0dfcad740d
# ╟─84a596c1-3a0a-4edd-a7ea-f3fcd8839c2a
# ╟─91b2718a-a921-4bae-99c0-91ec9c7e6479
# ╠═5a0b156e-b0b2-43b8-b15e-bcdb65c154bf
# ╟─a3270d8b-f167-4fb9-bcc5-c62c91b39649
# ╟─6528501f-3bb2-48a7-b24e-d32400033538
# ╟─0d0d6d2f-56b2-4e01-b9a4-6e469eac90ad
# ╟─76cfd5c8-6298-4ae2-9d2d-9e60deb12f49
# ╟─0bf40552-b627-4fa7-b3ba-1d8c7ed4a568
