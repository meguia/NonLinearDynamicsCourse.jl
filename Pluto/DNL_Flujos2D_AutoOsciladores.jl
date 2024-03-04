### A Pluto.jl notebook ###
# v0.19.39

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

# ╔═╡ 78c1fc6a-a145-4e30-9a9f-d4c66e1a591a
import Pkg; Pkg.add("ForwardDiff"), Pkg.add("Plots"), Pkg.add("PlutoUI"), Pkg.add("DifferentialEquations"), Pkg.add("IntervalRootFinding"), Pkg.add("StaticArrays"), Pkg.add("LinearAlgebra"), Pkg.add("Distances")

# ╔═╡ 17ccb720-4b26-4279-8a40-0e34a33d835d
push!(LOAD_PATH, joinpath(dirname(pwd()) , "src"))

# ╔═╡ 72fa8994-47dc-4e0c-80fd-021f5fc629e9
using NonLinearDynamicsCourse, PlutoUI, DifferentialEquations, Plots

# ╔═╡ 87c571dc-db6a-4cfe-9581-2b12d3f2ee91
md"""
# Oscilador de Van der Pol
 
Partimos de vuelta de las ecuaciones del oscilador armonico simple

$\dot{x} = y$

$\dot{y} = -\mu y - kx$
 
Recordamos que $y$ representa la velocidad del oscilador, $\dot{y}$ la aceleracion que es igual a la fuerza aplicada (suponemos masa igual a 1) que es lo que esta en el miembro derecho de la segunda ecuacion. En ese miembro $-kx$ representa la fuerza lineal elastica, mientras que $-\mu y$ es la friccion, una fuerza que siempre se opone a la velocidad ($\mu>0$) y que termina frenando el oscilador. En toda esta parte vamos a estudiar diferentes formas generales de esta **friccion** que de forma general vamos a ponerla como funcion de la posicion $x$ y la velocidad $y$:

$\dot{x} = y$

$\dot{y} = -F(x,y)y - kx$

Friccion lineal:

$F(x,y)=\mu$

Es interesante ver que $\mu$ se puede ver tambien como una resistencia.
De hecho las ecuaciones anteriores describen tambien otro sistema fisico muy estudiado, el circuito RLC:

<div>
<img src="../files/RLC.jpg" width="500px">
</div>

Para el caso en serie $\mu = R/L$ $k=1/LC$ y la variable $y$ corresponde a la corriente que circula por el circuito ($x$ vendria a ser la carga del capacitor).

Nos podemos preguntar ahora que pasaria si usaramos una "resistencia negativa" para valores pequeños de corriente (en realidad estrictamente de carga). Lo que queremos es que el termino $-\mu y$ se invierta para valores de $x$ pequeños. 

Por que? Porque de esta forma podemos evitar que las oscilaciones "mueran". Si la amplitud de la oscilacion $x$ se hace muy chica (el sistema se frena por la friccion o la resistencia) aparece una fuerza que va **a favor** de la velocidad (o de la circulacion de corriente) inyectandole energia al sistema, sin embargo para amplitudes de oscilacion grandes gana la disipacion (resistencia) que frena el sistema. 

De esta forma se alcanza un equilibrio en el que se producen auto-oscilaciones que no se extinguen. Estas oscilaciones en el espacio de fases se conocen como **ciclos limite** y son conjuntos invariantes (atractores o repulsores) como los puntos fijos. Notar la diferencia con las oscilaciones y orbitas concentricas del oscilador armonico sin friccion. A diferencia de este ultimo las oscilaciones de relajacion (o los ciclos limite en general) son atractoras, es decir que cualquier condicion cercana termina convergiendo a ellas. 


Como podemos escribir esta inversion de la friccion para valores de $x$ cercanos a cero? 
La forma mas simple es usar una nolinealidad cuadratica. Una parabola hundida en el eje tiene valores negativos cerca del $x=0$ y positivos para valores grandes de $x$ (tanto positivos como negativos). Por lo tanto podemos reemplazar a la friccion lineal por $\mu(x^2-1)$ . La resistencia (o friccion) en funcion de la posicion del oscilador (o de la carga en el capacitor) seria:

Friccion Van der Pol

$F(x,y) = \mu(x^2-1)$

Y el sistema de ecuaciones diferenciales queda escrito:

$\dot{x} = y$

$\dot{y} = \mu (1 -x^2)y - x$

La parte 'negativa' de la resistencia es $+\mu y$, una fuerza que va **a favor** de la velocidad. Estamos asumiendo que $\mu>0$. Mientras que queda un termino $- \mu x^2 y$ que siempre se va a oponer a la velocidad y gana para amplitudes grandes.

La idea de una 'resistencia negativa' no es caprichosa, algunos elementos electronicos de base de semiconductor y valvulares se comportan de esa forma (son activos). De hecho Van der Pol dedujo sus ecuaciones a partir del estudio de un circuito amplificador con un elemento valvular (triodo) en 1927. El mismo efecto se puede observar en tubos de neon y otros elementos con descarga en gases.


<div>
<img src="../files/Neon1.PNG" width="300px">
</div>

Una version mas elaborada de este sistema dio lugar al desarrollo del famoso VODER en 1939:

<div>
<img src="../files/Voder.PNG" width="500px">
</div>
"""

# ╔═╡ f19f478a-8428-4b25-9a10-9c4b85a3b096
md"""
## Puntos fijos y estabilidad
Vamos a calcular los puntos fijos primero por el metodo de las nulclinas.

La primera nuclina es trivial, una recta horizontal $y=0$ con lo cual todos los puntos fijos van a estar sobre el eje $x$

La segunda nulclina es mas complicada escrita como $y$ en funcion de $x$:

$y = \Large\frac{1}{\mu}\frac{x}{1-x^2}$

pero el unico punto que esta en las dos curvas es $(0,0)$ que es el unico punto fijo. 

Calculamos la matriz jacobiana


$
\begin{pmatrix}
0 & 1\\
-1-2\mu xy & \mu(1-x^2)
\end{pmatrix}
$

que evaluada en el unico punto fijo da

$
\begin{pmatrix}
0 & 1\\
-1 & \mu
\end{pmatrix}
$

Es facil calcular $Tr=\mu$ y $\Delta=1$. El determinante es siempre positivo y la traza tiene el signo del parámetro $\mu$ uqe asumimos positivo

Cuando $\mu>0$ el origen se hace repulsor, pero globalmente el flujo es atractor, para valores grandes de $x$ y de $y$ la disipacion no lineal $-x^2y$ domina y actua como un atractor global. Por lo tanto se va a formar un **ciclo limite estable** en el medio. Veamos como son las nulclinas y el flujo.
"""

# ╔═╡ cdc8339f-138f-4eb7-a88d-7ca495cc453f
#definimos la Ed para el oscilador de VanderPol
function vdp!(du,u,p,t)
    du[1] = u[2]
    du[2] = p[1]*(1.0-u[1]*u[1])*u[2]-u[1]
    du
end    

# ╔═╡ 0ea1cdee-a4a3-49b3-a599-7d3f3449855b
u0_arr=[[2.0;3.0],[-2.0;-3.0],[0.1;-0.1],[3.0;-2.0]]

# ╔═╡ 2f7d77f0-2eff-43ff-93d4-cba683e6030c
flux2d_nullclines(vdp!,u0_arr,10.0,[0.8];xlims=[-3,3],ylims=[-3,3],title="van der Pol")

# ╔═╡ 03ff660e-4b24-4721-b9f3-13d107600d4a
md"""
x0 $(@bind x0_vdp Slider(-1.0:0.1:2.0,default=0.1;show_value=true)) \
y0 $(@bind y0_vdp Slider(-1.0:0.1:1.0,default=0.1;show_value=true)) \
μ : $(@bind μ_vdp Slider(-0.1:0.01:10.0,default=-0.1;show_value=true)) \
tmax : $(@bind tmax_vdp Slider(10:10:300.0,default=10.0;show_value=true)) 
"""

# ╔═╡ 83131186-2dbb-4e2c-9d7a-6b2d829e2222
begin
	sol=solve(ODEProblem(vdp!,[x0_vdp,y0_vdp],(0,tmax_vdp),[μ_vdp]))
    p1=plot(sol,idxs=(0,1),xlabel="t",ylabel="x")
    p2=plot(sol,idxs=(1,2),legend=false,xlabel="x",ylabel="y")
    scatter!(p2,[sol.u[1][1]],[sol.u[2][1]])
    plot(p1,p2,layout=(1,2),size=(700,300),fmt=:png,title="van der Pol")
end	

# ╔═╡ 818e26ce-d164-4bad-8e31-7a1bf6be9ef3
md"""
Con la funcion `flux2d_animated` podemos hacer un gif animado de la evolucion del sistema. Recibe como argumentos:
- la funcion del campo vector
- el vector de parametros 
- la cantidad total de frames
- el paso temporal entre frames

Y como argumentos opcionales:
- Ngrid la cantidad de puntos de la grilla de condiciones iniciales (default Ngrid=10)
- fps (default fps=15) la cantidad de frames por segundo del gif
- xlmis,ylims los limites del grafico
- nullclines (default false) si agrega las nulclinas al grafico.
- fname (default empty) nombre del gif, por defecto le da el nombre de la funcion y el fps
"""

# ╔═╡ e3e59c04-2817-44be-8df0-3e42f77cda58
flux2d_animated(vdp!,[1.0],100,0.02;Ngrid=10,xlims=[-3,3],ylims=[-15,15])

# ╔═╡ 136021df-bc6b-4e0c-b136-51ecbe46c2ff
md"""
## Osciladores de relajacion. Transformacion de Lienard

A diferencia de los puntos fijos no hay un método general para establecer la existencia de ciclos límites, es decir de órbitas cerradas **aisladas** que funcionan como conjunto atractor (o repulsor). Existen criterios que por lo general son difíciles de aplicar y muy frecuentemente la existencia de ciclos se puede presuponer a partir del análisis del flujo global o a través de ciertas transformaciones de coordenadas. 

Vamos a ver un ejemplo para el oscilador de van der Pol para el caso en el que el parámetro $\mu$ tiene un valor alto. En ese caso el sistema se aproxima al comportamiento de un **oscilador de relajación**. Los osciladores de relajación surgen de la alternancia entre un proceso lento de carga y uno rápido de descarga (por ejemplo los picos en el comportamiento eléctrico de las neuronas) y suceden en dos escalas temporales diferentes. Si se puede hacer un cambio de coordenadas de forma tal que esas dos escalas de tiempo estén separadas en las nuevas variables (es decir que tengamos una variable "rápida" y otra "lenta") el ciclo límite va a formarse por una trayectoria que salta (o relaja) en la variable rápida de una rama a otra de la nulclina sobre la que se arrastra siguiendo la ecuación de la variable lenta.

En el caso del oscilador de van der Pol se puede hacer un cambio de coordenadas no lineal (sin meternos en los detalles ) que transforma la forma anterior a una oscilador de relajación con variables rápida $(x_1)$ y lenta $(x_2)$. El cambio es el siguiente:

$x_1 = \mu^{-1/2}x$

$y_1 = \mu^{-3/2}(\mu x - x^3/3 -y)$

En las nuevas variables el sistema queda escrito como 

$\dot{x_1} = \mu(x_1-x_1^3/3-y_1)$

$\dot{y_1} = x_1/\mu$

como antes el unico punto fijo esta en $(0,0)$ pero la primera nulclina se puede escribir como una cubica que tiene siempre forma de "N":

$y_1 = x_1 - x_1^2/3$

y la otra es una recta vertical en $x_1=0$

Pero lo interesante pasa cuando $\mu$ es grande. Como esta dividiendo a la variacion de $y_1$ y multiplicando a la variacion de $x_1$, podemos imaginar que el flujo va a ser mucho mas rapido en la direccion horizontal $\dot{x_1}>>\dot{x_2}$, salvo cuando se aproxima a la nulclina en forma de "N" (porque ahi se hace $\dot{x_1}=1$). 
Entonces tenemos un flujo rapido horizontal que va a parar a la nulclina en forma de "N". 

Toda la zona que esta arriba de la N ($y_1>x_1-x_1^3/3$ hace $\dot{x_1}$ sea negativo) fluye hacia la izquierda y lo que esta abajo de la N ($y_1<x_1-x_1^3/3$ hace $\dot{x_1}$ sea positivo) fluye a la derecha.

Esto hace que en terminos del flujo horizontal la parte del medio de la N sea inestable (el flujo se aleja de ella) y las ramas de los costados sean estables (el flujo es atraido hacia ellas)

Una vez que la trayectoria esta cerca de la nulclina N va a fluir lentamente siguiendo el signo de $x_1$: hacia abajo a la izquierda del eje vertical y hacia arriba a la derecha del eje. 

Cuando la orbita llega al "codo" de la N se encuentra con la rama inestable de la nulclina y 'salta' (relaja) a la otra rama estable. Esto es lo que se conoce como un oscilador de relajacion.

Atencion porque esta forma sirve para ver mejor el ciclo limite pero por la transformacion de coorenadas la bifurcacion de Hopf es 'degenerada' ($\mu=0$ implicaria ademas que $\dot{y_1}$ se hiciese infinito).
"""

# ╔═╡ 13cfddda-608e-4cec-8acb-2f129c3cee93
#definimos la Ed para el oscilador de VanderPol en la version de Lienard (oscilador de relajacion)
function lienard!(du,u,p,t)
    du[1] = p[1]*(u[1]*(1.0-u[1]*u[1]/3.0)-u[2])
    du[2] = u[1]/p[1]
end    

# ╔═╡ d1c1644b-17e1-4c1d-aa2b-60cae1ebf5cd
flux2d_nullclines(lienard!,[2.0;0.1],20.0,[3.0],xlims=[-3,3],ylims=[-2,2],title="Oscilador de Relajacion")

# ╔═╡ 465c2e99-b477-47a8-8503-4f356bf1dd7b
md"""
# Auto-Oscilador de Rayleigh. Lengueta. 

El oscilador de van der Pol lo habiamos presentado como el auto-oscilador con 'friccion negativa' mas simple posible porque la no linealidad era cuadratica y estaba solo en la variable $x$. Esa dependencia de la resistencia en $x$ aparecia de forma natural en circuitos valulares (tristores) y lo presentamos por motivos historicos.

Pero la no linealidad que da una zona de friccion negativa puede ser una funcion de la velocidad $y$ en lugar de $x$. De hecho esta forma aparece naturalmente cuando la fuerza que da impulso al oscilador cuando tiene poca amplitud es de origen aerodinamico. Entonces en su forma mas simple, con una no linealidad cuadratica podemos proponer la friccion no lineal.

Friccion Rayleigh:

$C(x,y)=y^2-\mu$

Comparar con la friccion de Van der Pol del principio. Esta forma la llamamos de Rayleigh porque es quien en su Teoria del Sonido en 1865 propone este modelo para la oscilacion de una lengueta de clarinete:

$\dot{x} = y$

$\dot{y} = -(y^2-\mu)y-kx$

donde $k$ determina la frecuencia de oscilacion. El sistema anterior se puede escribir entonces, distribuyendo los dos terminos de la fraccion:

$\dot{x} = y$

$\dot{y} = \mu y-y^3-kx$

De forma analoga que en el modelo de Van der Pol el termino $y^3$ representa la friccion no lineal que atenua las oscilaciones para amplitudes grandes y $\mu y$ es la 'friccion negativa' que actua para valores pequeños de $y$ como fuerza restitutiva y representa la accion del clarinetista generando una inestabilidad con el flujo del aire contra la lengueta.

El calculo de las nulclinas los puntos fijos y la estabilidad de los mismo queda para la practica. Pero se puede notar al jugar con los parametros que surge un ciclo límite estable al igual que en el Van der Pol. Para que valor de $\mu$?
"""

# ╔═╡ b141b69f-2c2c-4928-94fd-3f1ae92d38ac
# Escribimos la ecuacion de la lengueta
function reed!(du,u,p,t)
    (μ,k) = p
    du[1] = u[2]
    du[2] = u[2]*(μ-u[2]*u[2])-k*u[1]
    du
end    

# ╔═╡ 835a9bda-3708-41af-bdab-1b8e72aac4aa
md"""
x0 $(@bind x0_reed Slider(-1.0:0.1:2.0,default=0.1;show_value=true)) \
y0 $(@bind y0_reed Slider(-1.0:0.1:1.0,default=0.1;show_value=true)) \
μ : $(@bind μ_reed Slider(-0.1:0.01:10.0,default=-0.1;show_value=true)) \
k : $(@bind k_reed Slider(0.1:0.01:3.0,default=-0.1;show_value=true)) \
tmax : $(@bind tmax_reed Slider(10:10:300.0,default=10.0;show_value=true)) 
"""

# ╔═╡ df7cab14-abe1-4f93-aec8-5422fea44fde
flux2d_nullclines(reed!,[x0_reed;y0_reed],tmax_reed,[μ_reed , k_reed],xlims=[-3,3],ylims=[-2,2],title="Lengüeta (Rayleigh)")

# ╔═╡ 163d352d-ce42-4828-a418-520033c23719
# Escribimos la ecuacion de la lengueta
function vreed!(du,u,p,t)
    (μ,k,v0) = p
	v = u[2] - v0
    du[1] = u[2]
    du[2] = v*(μ-v*v)-k*u[1]
    du
end    

# ╔═╡ e8947aad-3033-45e1-b9be-c85950b82fef
md"""
x0 $(@bind x0_vreed Slider(-1.0:0.1:2.0,default=0.1;show_value=true)) \
y0 $(@bind y0_vreed Slider(-1.0:0.1:1.0,default=0.1;show_value=true)) \
μ : $(@bind μ_vreed Slider(-0.1:0.01:10.0,default=-0.1;show_value=true)) \
v0 : $(@bind v0_vreed Slider(-1.0:0.01:1.0,default=-0.1;show_value=true)) \
k : $(@bind k_vreed Slider(0.1:0.01:3.0,default=-0.1;show_value=true)) \
tmax : $(@bind tmax_vreed Slider(10:10:300.0,default=10.0;show_value=true)) 
"""

# ╔═╡ ea88f136-6d7d-4b0a-b4d9-2f4e898930af
flux2d_nullclines(vreed!,[x0_vreed;y0_vreed],tmax_vreed,[μ_vreed , k_vreed, v0_vreed],xlims=[-3,3],ylims=[-2,2],title="Lengüeta (Rayleigh) Modificada")

# ╔═╡ 5a415de4-126d-462c-9eb0-ee06500dae22
md"""
# Oscilador Frotado (violin puntual)

Otro sistema con auto-oscilaciones simples fue propuesto (tambien por Rayleigh en 1877!) para modelar la accion de **slip & stick** del arco contra la cuerda del violin, pero se puede aplicar a un monton de sistemas que generan autooscilaciones a partir de la friccion.

<div>
<img src="../files/Conveyor.PNG" width="300px">
</div>


El modelo propuesto era similar al de la figura. Una masa unida a un resorte esta apoyada sobre una cinta transportadora con friccion que se mueve con velocidad constante hacia la derecga. Al principio el rozamiento estatico hace que la masa este adherida (momento **stick**) a la cinta y ejerce una fuerza que iguala a la del resorte. Pero la friccion estatica tiene un valor maximo y cuando el resorte esta muy estirado no puede superar la fuerza elastica y la masa es arrastrada por el resorte y desliza con rozamiento dinamico sobre la cinta hacia la izquierda (momento **slip**) y puede llegar por inercia incluso a comprimir un poco el resorte hasta que la masa se frena y queda enganchada de vuelta por el rozamiento estatico y se repite el proceso.

Si bien se trata de un sistema simple, la forma funcional de la friccion (que tiene que ser funcion de la diferencia de velocidad entre la masa y la cinta es decir si desliza o no) no puede ser algo tan simple como una cuadratica o una cubica porque tiene que cambiar de signo bruscamente, ya que la friccion tiene que ser maxima para deslizamientos bajos e ir decreciendo para deslizamientos mas rapidos. La forma clasica es algo asi:

<div>
<img src="../files/Friction.PNG" width="300px">
</div>

donde $\dot{x}-v$ es el 'deslizamiento', es decir la diferencia de velocidades entre la masa y la cinta. Cuando la masa esta adherida a la cinta $C$ puede tomar todos los valores en la vertical hasta un valor maximo hacia un lado y hacia el otro y luego 'salta' al rozamiento dinamico que es menor a medida que el deslizamiento es mas cada vez rapido.

La forma funcional que esta representada arriba para el desplazamiento $d=\dot{x}-v$  es:

$C(d)=sign(d) e^{-2|d|}$

Ese salto brusco de la funcion signo trae problemas numericos, lo regularizamos con la funcion arco tangente del desplazamiento dividido por un numero pequeño. el arcotangente es un escalon mas suave que el signo (si el valor se hace muy pequeño se va haciendo cada vez mas parecido a la funcion signo).

$C_{bow}(d)=arctan(d/\epsilon) e^{-2|d|}$

Con todos estos elementos la friccion del Arco propuesta por Rayleigh y el sistema final queda escrito

$\dot{x}=y$

$\dot{y}=-\mu C_{bow}(y-v) -x$
"""

# ╔═╡ a86dcb92-ffc6-44f5-9be0-d787b592cb44
friction(x) = atan(x/0.05)*exp(-2*abs(x))

# ╔═╡ 678cc7eb-8e86-47ec-8176-35d30eff5587
function bow!(du,u,p,t)
    du[1]=u[2]
    du[2]=-p[1]*friction(u[2]-p[2])-u[1]
    du
end    

# ╔═╡ 66497e63-8efd-472b-aa35-50c826917eb2
md"""
x0 $(@bind x0_bow Slider(-1.0:0.1:2.0,default=0.1;show_value=true)) \
y0 $(@bind y0_bow Slider(-1.0:0.1:1.0,default=0.1;show_value=true)) \
μ : $(@bind μ_bow Slider(0.01:0.01:1.0,default=0.1;show_value=true)) \
v0 : $(@bind v0_bow Slider(-1.0:0.01:1.0,default=-0.1;show_value=true)) \
tmax : $(@bind tmax_bow Slider(10:10:300.0,default=10.0;show_value=true)) 
"""

# ╔═╡ 35f0087d-28fc-42b8-baeb-c18df8a149fd
flux2d_nullclines(bow!,[x0_bow;y0_bow],tmax_bow,[μ_bow,v0_bow];xlims=[-3,3],ylims=[-2,2],title="Cuerda Frotada")

# ╔═╡ fdb0d8bc-d780-4b81-9e54-19510d70ed76
begin
	sol2 = solve(ODEProblem(bow!, [x0_bow;y0_bow],tmax_bow,[μ_bow,v0_bow]));
    p1b = plot(sol2)
    p2b = plot(sol2,idxs=(1,2),arrow=true)
    plot(p1b,p2b,layout=(1,2),size = (900,450),title="Cuerda Frotada")
end	

# ╔═╡ 7fb98f93-be03-43d3-895c-c72f58f09d27
md"""
# Oscilador de Duffing van der Pol

En todos los osciladores anteriores la fuerza de restitucion era lineal, pero podemos considerar osciladores de forma mas general como:

$\dot{x} = y$

$\dot{y} = -C(x,y)y - K(x)$

en el caso del oscilador armonico $C(x,y)=c$ y $K(x)=kx$. En el oscilador de Van der Pol y en de Rayleigh cambia la forma de $C(x,y)$ pero la restitucion sigue siendo lineal. 

El oscilador de Duffing propone una fuerza de restitucion cubica: $K(x)=x^3-\beta x$.

Podemos combinar la restitucion cubica de Duffing con la friccion negativa de Van der Pol para tener un auto-oscilador no lineal con mas variedad de comportamiento. Tenemos entondes el modelo de Duffing-Van der Pol:

$\dot{x} = y$

$\dot{y} = \mu y -x^2y +\beta x - x^3$

"""

# ╔═╡ 216ef632-f17a-4f28-9aad-04364803bcb7
#definimos la Ed para el oscilador de Duffing-VanderPol
function duffing_vanderpol!(du,u,p,t)
    du[1] = u[2]
    du[2] = (p[1]-u[1]*u[1])*u[2]+u[1]*(p[2]-u[1]*u[1])
    du
end    

# ╔═╡ b81e8195-10ac-467f-b99f-c881a4f3d681
md"""
x0 $(@bind x0_dvdp Slider(-1.0:0.1:2.0,default=0.1;show_value=true)) \
y0 $(@bind y0_dvdp Slider(-1.0:0.1:1.0,default=0.1;show_value=true)) \
μ : $(@bind μ_dvdp Slider(0.01:0.01:1.0,default=0.1;show_value=true)) \
β : $(@bind β_dvdp Slider(-1.0:0.01:1.0,default=-0.1;show_value=true)) \
tmax : $(@bind tmax_dvdp Slider(10:10:300.0,default=10.0;show_value=true)) 
"""

# ╔═╡ fe8fe3ec-34f2-40e8-a48e-c5d9cbc50985
flux2d_nullclines(duffing_vanderpol!,[x0_dvdp;y0_dvdp],tmax_dvdp,[μ_dvdp,β_dvdp]; npts=41,xlims=[-3,3],ylims=[-2,2],title="Duffing van der Pol")

# ╔═╡ 362650ef-1899-4e53-93c0-c5761bdf47b8
begin
	sol3 = solve(ODEProblem(duffing_vanderpol!, [x0_dvdp;y0_dvdp],tmax_dvdp,[μ_dvdp,β_dvdp]));
    p1c = plot(sol3)
    p2c = plot(sol3,idxs=(1,2),arrow=true)
    plot(p1c,p2c,layout=(1,2),size = (900,450),title="Duffing Van der Pol")
end	

# ╔═╡ 4f04d9a1-3e58-4ead-b5c4-daf4d88e9f29
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
# ╠═17ccb720-4b26-4279-8a40-0e34a33d835d
# ╠═78c1fc6a-a145-4e30-9a9f-d4c66e1a591a
# ╠═72fa8994-47dc-4e0c-80fd-021f5fc629e9
# ╟─87c571dc-db6a-4cfe-9581-2b12d3f2ee91
# ╟─f19f478a-8428-4b25-9a10-9c4b85a3b096
# ╠═cdc8339f-138f-4eb7-a88d-7ca495cc453f
# ╠═0ea1cdee-a4a3-49b3-a599-7d3f3449855b
# ╠═2f7d77f0-2eff-43ff-93d4-cba683e6030c
# ╠═03ff660e-4b24-4721-b9f3-13d107600d4a
# ╠═83131186-2dbb-4e2c-9d7a-6b2d829e2222
# ╟─818e26ce-d164-4bad-8e31-7a1bf6be9ef3
# ╠═e3e59c04-2817-44be-8df0-3e42f77cda58
# ╟─136021df-bc6b-4e0c-b136-51ecbe46c2ff
# ╠═13cfddda-608e-4cec-8acb-2f129c3cee93
# ╠═d1c1644b-17e1-4c1d-aa2b-60cae1ebf5cd
# ╟─465c2e99-b477-47a8-8503-4f356bf1dd7b
# ╠═b141b69f-2c2c-4928-94fd-3f1ae92d38ac
# ╠═835a9bda-3708-41af-bdab-1b8e72aac4aa
# ╠═df7cab14-abe1-4f93-aec8-5422fea44fde
# ╠═163d352d-ce42-4828-a418-520033c23719
# ╠═e8947aad-3033-45e1-b9be-c85950b82fef
# ╠═ea88f136-6d7d-4b0a-b4d9-2f4e898930af
# ╟─5a415de4-126d-462c-9eb0-ee06500dae22
# ╠═a86dcb92-ffc6-44f5-9be0-d787b592cb44
# ╠═678cc7eb-8e86-47ec-8176-35d30eff5587
# ╟─66497e63-8efd-472b-aa35-50c826917eb2
# ╠═35f0087d-28fc-42b8-baeb-c18df8a149fd
# ╠═fdb0d8bc-d780-4b81-9e54-19510d70ed76
# ╟─7fb98f93-be03-43d3-895c-c72f58f09d27
# ╠═216ef632-f17a-4f28-9aad-04364803bcb7
# ╟─b81e8195-10ac-467f-b99f-c881a4f3d681
# ╠═fe8fe3ec-34f2-40e8-a48e-c5d9cbc50985
# ╠═362650ef-1899-4e53-93c0-c5761bdf47b8
# ╟─4f04d9a1-3e58-4ead-b5c4-daf4d88e9f29
