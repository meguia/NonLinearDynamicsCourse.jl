### A Pluto.jl notebook ###
# v0.19.40

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

# ╔═╡ 1d36f38e-d842-47f1-8577-a16bb9a73d76
using Pkg; Pkg.activate(joinpath(dirname(pwd()) , "src"))

# ╔═╡ bc0b7062-f045-11ee-260f-f99d879c90be
push!(LOAD_PATH, joinpath(dirname(pwd()) , "src"))

# ╔═╡ 69f343ec-574b-4252-9134-b30bbd74effb
using NonLinearDynamicsCourse, PlutoUI, DifferentialEquations, Plots

# ╔═╡ 787ce7fe-8a11-48c3-aca9-9a4c0ac4dbd2
md""" 
# Bautin (Hopf Generalizada)

La forma normal de la bifurcacion de Bautin (codimension 2) es:

$\dot{\mathtt{x}} = \beta_1 \mathtt{x}- \omega \mathtt{y} + \beta_2 \mathtt{x} (\mathtt{x}^2+\mathtt{y}^2) + σ\mathtt{x} (\mathtt{x}^2+\mathtt{y}^2)^2$

$\dot{\mathtt{y}} = \omega \mathtt{x} + \beta_1 \mathtt{y} + \beta_2 \mathtt{y} (\mathtt{x}^2+\mathtt{y}^2) + σ\mathtt{y} (\mathtt{x}^2+\mathtt{y}^2)^2$

![Bautin](http://www.scholarpedia.org/w/images/d/d0/Bautin.gif)

"""

# ╔═╡ c2177a4c-cff4-4c7a-86f4-629b60e13a63
md"""
## Bautin de un sistema estandar

Planteamos un sistema estandar siguiendo la estrategia de la reduccion de la forma normal de Hopf para elegir los terminos no lineales.

$\dot{x} = y$

$\dot{y} = -kx + 2\mu y + \sigma_{112} x^2y + \sigma_{222} y^3 + \xi O(5)$


### Autovalores y diagonalizacion de la parte Lineal

Escribimos el sistema anterior como 

$\dot{\bf{x}}=A\bf{x}+\bf{F}(\bf{x})$

con la parte lineal:

$A=\begin{pmatrix}
0 & 1\\
-k & 2\mu
\end{pmatrix}$

cuyos autovalores son los complejos conjugados: $\lambda_1 = \lambda$ y $\lambda_1 = \bar{\lambda}$ con:

$\lambda = \mu + i \sqrt{k-\mu^2}$

con autovectores $q_1=(1,\lambda)^T$ , $q_2=(1,\bar{\lambda})^T$. 

El procedimiento ahora es introducir una variable compleja $z$ usando la transpuesta de A (aunque tambien se puede proyectar sobre la base real de A). La matriz transpuesta:

$A^T = \begin{pmatrix}
0 & -k\\
1 & 2\mu
\end{pmatrix}$

y sus autovectores: $p_1=(-\bar{\lambda},1)^T$ y $p_2=(-\lambda,1)^T$. 

Podemos verificar que los autovalores son ortogonales (recordando que el producto escala complejo conjuga a derecha)

$<p_1,q_1>=-\bar{\lambda}+\bar{\lambda}=0$

y podemos normalizar $p_1$ de forma tal que $<p_1,q_2>=-\bar{\lambda}+\lambda$ sea igual a uno.

El autovector de la transpuesta normalizado es entonces

$\hat{p_1} = \frac{1}{\lambda-\bar{\lambda}} \begin{pmatrix}-\bar{\lambda}\\1\end{pmatrix}= \frac{1}{i2\omega}\begin{pmatrix}-\bar{\lambda}\\1\end{pmatrix}$

y la nueva variable compleja definida como el producto escalar $z = <\hat{p_1},\bf{x}>$:

$z = \frac{y-\bar{\lambda}x}{\lambda-\bar{\lambda}} = \frac{1}{i2\omega}(y-\bar{\lambda}x)$

donde definimos $\omega=\sqrt{k-\mu^2}$ 

podemos tener las transformaciones inversas como 

$x = z + \bar{z} = 2Re(z)$

$y = \lambda z + \bar{\lambda}\bar{z} = 2\mu Re(z) - 2\omega Im(z)$



se puede verificar que

$\dot{z} = \frac{1}{i2\omega}(\dot{y} - \bar{\lambda}\dot{x}) = \frac{1}{i2\omega}(-kx + 2\mu y - \mu y + i \omega y + \sigma O(3) + \xi O(5)) = \frac{1}{i2\omega}(\lambda y - \lambda\bar{\lambda}x + \sigma O(3) + \xi O(5))=\lambda z + \frac{1}{i2\omega}(\sigma O(3) + \xi O(5))$

donde usamos el hecho de que $\lambda \bar{\lambda} = k$

Si elegimos como nuevas variables la parte real e imaginaria de $z$ la parte lineal nos queda con la matriz:

$\begin{pmatrix}
\mu & -\omega\\
\omega & \mu
\end{pmatrix}$

Que es la parte lineal de la forma normal de la Hopf. 

Podemos evaluar las condiciones de genericidad  transversalidad y de no degeneracion y. La primera es sencilla porque el parametro original $\mu$ es identico a la parte real del autovalor diagonalizado.

Con respecto a la condicion de no degeneracion (que es la que se tiene que violar en la Bautin), tenemos que evaluar el primer coeficiente de Lyapunov

## Genericidad de la Hopf y Terminos no lineales

La estrategia usual es partir del desarrollo de Taylor multilineal en la bifurcacion:

$F(\bf{x}) = \frac{1}{2} B(\bf{x},\bf{x}) + \frac{1}{6} C(\bf{x},\bf{x},\bf{x}) + O(4)$

y decidir que terminos vamos a incluir en funcion de como se incluyan en el coeficiente de Lyapunov dependiendo de la proyeccion dada por los autovectores $q$ y $p$ de la parte lineal. En nuestro caso no vamos a incluir terminos cuadraticos y solo terminos cubicos impares en $y$, 
es decir 

$F(\bf{x}) = \begin{pmatrix}
0\\
\sigma_{112} x_1^2x_2 + \sigma_{222} x_2^3
\end{pmatrix}$

donde para mayor consistencia con lo que viene usamos la notacion $\bf{x}=(x_1,x_2)$ en lugar de las variables $x$ $y$
Pero para poder hacer la proyeccion tenemos que escribir la funcion multilineal:

$C_{\alpha\beta\gamma} =  \begin{pmatrix}
0\\
2\sigma_{112} (\alpha_1\beta_1\gamma_2+\alpha_1\beta_2\gamma_1+\alpha_2\beta_1\gamma_1) + 6\sigma_{222} \alpha_2\beta_2\gamma_2
\end{pmatrix}$

El coeficiente de Lyapunov en la bifurcacion cuando no hay terminos cuadraticos queda expresado como

$l_1 = \frac{1}{2\omega}Re(<\hat{p}_1,C(q,q,\bar{q})>)$

Evaluamos la segunda componente de $C$ en los autovectores $q,q,\bar{q}$ recordando que $q=(1,\lambda)^T$ 

$C_2(q,q,\bar{q}) = 2\sigma_{112} (2\lambda+\bar{\lambda}) + 6 \sigma_{222} (\lambda^2\bar{\lambda})= 2\sigma_{112} i\omega + 6 \sigma_{222} i \omega^3$

Recordar que estamos evaluando en la bifuracacion $\mu=0$ entonces $\lambda = - \bar{\lambda} = i\omega$.

El ultimo paso es proyectar con $p_1$ y tomar la parte real (que es la unica que sobrevive por la normalizacion)

$<\hat{p}_1,C(q,q,\bar{q})> = \left< \frac{1}{i2\omega}\begin{pmatrix}-\bar{\lambda}\\1\end{pmatrix} ,\begin{pmatrix} 0\\
2i\omega(\sigma_{112}+3\sigma_{222}\omega^2)\end{pmatrix}\right> = \sigma_{112}+3\sigma_{222}\omega^2$

si bien hay una normalizacion adicional, nos interesa el signo de $l_1$ y este esta dado por el signo de:

$\sigma_{112}+3\sigma_{222}k$

En el ejemplo de Van der Pol de quinto orden ambos terminos estan presentes y con el mismo valor, lo cual hace que la zona biestable dependa del parametro $k$. Para independizarnos de esto y para que quede el sistema lo mas sencillo posible nos quedamos solo con el termino $\sigma_{112}$ y pasamos a llamar simplemente $\sigma$ a ese parametro. 


Queda por determinar entonces cual es el termino (o los terminos) de quinto orden que sobreviven. En este caso es mas dificil y elegimos tomar el de la forma normal con la potencia cuarta del modulo de z ya que es el que va a garantizar la estabilidad global siempre que este con signo negativo. El sistema estandar final por lo tanto queda escrito (con $\mu'=2\mu$ del parametro original pero luego vamos a quitar la pima) como:

$\dot{x} = y$

$\dot{y} = -kx + \mu' y + \sigma x^2y  + \xi y (x^2+y^2)^2$

"""

# ╔═╡ 0d7be9cf-763d-4781-b747-ee4f99fabeee
function bautin!(du,u,p,t)
	(μ,k,σ,ξ) = p
	du[1] = u[2]
	du[2] = -k*u[1] + u[2]*(μ + σ*u[1]^2 + ξ*(u[1]^2+u[2]^2)^2)
end

# ╔═╡ c8148366-9fab-4ef4-8ac8-0bb6421ce68e


# ╔═╡ ef53daf7-7d49-4763-8f3c-01aec8a105f1
md"""
# Bautin acopladas (Doble Hopf)

$\dot{x_1}  =  y_1$
$\dot{y_1}  =  -k_1x_1 + \mu_1 y_1 + \sigma x_1^2y_1  + c_{21} x_2^2y_1  + \xi y_1 (x_1^2+y_1^2)^2$
$\dot{x_2}  =  y_2$
$\dot{y_2}  =  -k_2x_2 + \mu_2 y_2 + \sigma x_2^2y_2  + c_{12} x_1^2y_2  + \xi y_2 (x_2^2+y_2^2)^2$


"""

# ╔═╡ bdcf36f8-a090-4a8a-a3a3-838defcc5b78
function dhopf!(du,u,p,t)
	(μ1,μ2,k1,k2,σ,c12,c21,ξ) = p
	du[1] = u[2]
	du[2] = -k1*u[1] + u[2]*(μ1 + σ*u[1]^2 + c21*u[3]^2 + ξ*(u[1]^2+u[2]^2)^2)
	du[3] = u[4]
	du[4] = -k2*u[3] + u[4]*(μ2 + σ*u[3]^2 + c12*u[1]^2 + ξ*(u[3]^2+u[4]^2)^2)
end

# ╔═╡ 53b50ba4-a968-4b6c-8496-031b5d72f558
html"""
<style>
main {
    max-width: 1000px;
}
input[type*="range"] {
	width: 40%;
}
</style>
"""

# ╔═╡ c59e2248-b18b-4ea1-a8ea-9b300de44f13
sp = html"&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp";

# ╔═╡ f3ac4ce9-d182-4b75-875e-274ee2a7f470
begin 
	md"""
	x0 $(@bind x0 Slider(-1.0:0.01:2.0,default=0.1;show_value=true)) $sp 
	y0 $(@bind y0 Slider(-1.0:0.1:1.0,default=0.1;show_value=true)) \
	μ : $(@bind μ Slider(-1.0:0.01:1.0,default=0.1;show_value=true)) $sp 
	k : $(@bind k Slider(0.01:0.01:2.0,default=1.0;show_value=true)) \
	σ : $(@bind σ Slider(-1.0:0.01:1.0,default=0.1;show_value=true)) $sp 
	ξ : $(@bind ξ Slider(-1.0:0.01:-0.01,default=-0.1;show_value=true)) \
	tmax : $(@bind tmax Slider(10:10:3000.0,default=10.0;show_value=true)) 
	"""
end	

# ╔═╡ cc1e8cbc-1088-4ae0-bad0-8dcc83b98f22
begin
	sol = solve(ODEProblem(bautin!, [x0;y0],tmax,[μ,k,σ,ξ]));
	p1 = plot(sol)
    p2 = plot(sol,idxs=(1,2),arrow=true)
	plot(p1,p2,layout=(1,2),size = (900,450),title="Bautin")
end	

# ╔═╡ 5621e2f8-69b1-4bcc-a347-704d6e1bc235
flux2d_nullclines(bautin!,[x0;y0],tmax,[μ,k,σ,ξ],xlims=[-5,5],ylims=[-3,3],npts=60,title="Bautin")

# ╔═╡ 944a42df-4696-41f1-bc48-da1adaac122b
md"""
x10 $(@bind x10 Slider(-1.0:0.01:2.0,default=0.1;show_value=true)) $sp
x20 $(@bind x20 Slider(-1.0:0.1:1.0,default=0.1;show_value=true)) \
μ1 : $(@bind μ1 Slider(-1.0:0.01:1.0,default=-0.1;show_value=true)) $sp
μ2 : $(@bind μ2 Slider(-1.0:0.01:1.0,default=-0.1;show_value=true)) \
k1 : $(@bind k1 Slider(0.01:0.01:2.0,default=1.0;show_value=true)) $sp
k2 : $(@bind k2 Slider(0.01:0.01:2.0,default=1.0;show_value=true)) \
c12 : $(@bind c12 Slider(-1.0:0.01:1.0,default=0.0;show_value=true)) $sp
c21 : $(@bind c21 Slider(-1.0:0.01:1.0,default=0.0;show_value=true)) \
σ : $(@bind σ2 Slider(-1.0:0.01:1.0,default=-0.1;show_value=true)) $sp
tmax : $(@bind tmax2 Slider(10:10:3000.0,default=10.0;show_value=true)) 	
"""

# ╔═╡ b07e11bc-635f-4ab6-9ec3-b4c89d0347be
begin
	sol1 = solve(ODEProblem(dhopf!, [x10;0.0;x20;0.0],tmax2,[μ1,μ2,k1,k2,σ2,c12,c21,-0.1]));
	p1b = plot(sol1,idxs=(1,2),arrow=true)
	p2b = plot(sol1,idxs=(3,4),arrow=true)
	plot(p1b,p2b,layout=(1,2),size = (900,450),title="Double Hopf")
end	

# ╔═╡ 32cec679-0d93-4d26-a736-e0e6c1058f89
begin
	plot(sol1,idxs=(0,1),size=(1000,300),label="x1")
	plot!(sol1,idxs=(0,3),size=(1000,300),label="x2")
end	

# ╔═╡ eebdc7e9-1ecc-481f-936e-9ea0da2332ea
 
md"""
x0 $(@bind x0_reed Slider(-1.0:0.01:2.0,default=0.1;show_value=true)) $sp 
y0 $(@bind y0_reed Slider(-1.0:0.1:1.0,default=0.1;show_value=true)) \
μ : $(@bind μ_reed Slider(-1.0:0.01:1.0,default=0.1;show_value=true)) $sp 
k $(@bind k_reed Slider(0.0:0.1:1.0,default=0.1;show_value=true)) \
v0 : $(@bind v0 Slider(-1.0:0.01:1.0,default=0.1;show_value=true)) $sp
tmax : $(@bind tmax_reed Slider(10:10:300,default=0.1;show_value=true)) 
"""	

# ╔═╡ Cell order:
# ╠═bc0b7062-f045-11ee-260f-f99d879c90be
# ╠═1d36f38e-d842-47f1-8577-a16bb9a73d76
# ╠═69f343ec-574b-4252-9134-b30bbd74effb
# ╟─787ce7fe-8a11-48c3-aca9-9a4c0ac4dbd2
# ╟─c2177a4c-cff4-4c7a-86f4-629b60e13a63
# ╠═0d7be9cf-763d-4781-b747-ee4f99fabeee
# ╠═cc1e8cbc-1088-4ae0-bad0-8dcc83b98f22
# ╠═f3ac4ce9-d182-4b75-875e-274ee2a7f470
# ╠═5621e2f8-69b1-4bcc-a347-704d6e1bc235
# ╠═c8148366-9fab-4ef4-8ac8-0bb6421ce68e
# ╟─ef53daf7-7d49-4763-8f3c-01aec8a105f1
# ╠═bdcf36f8-a090-4a8a-a3a3-838defcc5b78
# ╟─b07e11bc-635f-4ab6-9ec3-b4c89d0347be
# ╟─944a42df-4696-41f1-bc48-da1adaac122b
# ╠═32cec679-0d93-4d26-a736-e0e6c1058f89
# ╟─eebdc7e9-1ecc-481f-936e-9ea0da2332ea
# ╟─53b50ba4-a968-4b6c-8496-031b5d72f558
# ╟─c59e2248-b18b-4ea1-a8ea-9b300de44f13
