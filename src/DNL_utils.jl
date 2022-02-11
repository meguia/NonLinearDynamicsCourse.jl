using Plots
using LinearAlgebra
using DifferentialEquations
using ForwardDiff
using IntervalRootFinding
using StaticArrays

function plot_nullclines(f,p;xlims=[-1.0,1.0],ylims=[-1.0,1.0],npts=30,regions=true)
    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    xgrid = xlims[1]:xrange/npts:xlims[2]
    ygrid = ylims[1]:xrange/npts:ylims[2]
    u = [[x;y] for y in ygrid, x in xgrid]
    Z = u*0
    f.(Z,u,(p,),0)
    p1 = plot(legend=:none)
    if regions
        contourf!(p1,xgrid,ygrid,getindex.(Z,1),levels=[-1000,0,1000],alpha=0.8, c=[:red,:green],label="")
        contourf!(p1,xgrid,ygrid,getindex.(Z,2),levels=[-1000,0,1000],alpha=0.4, c=[:blue,:yellow],label="")
        contour!(p1,xgrid,ygrid,getindex.(Z,1),levels=[0],c=:red,label="")
        contour!(p1,xgrid,ygrid,getindex.(Z,2),levels=[0],c=:blue,label="")
    else
        contour!(p1,xgrid,ygrid,getindex.(Z,1),levels=[0],c=:gray,alpha=0.5,label="")
        contour!(p1,xgrid,ygrid,getindex.(Z,2),levels=[0],c=:gray,alpha=0.5,label="")
    end
    p1
end

function solve_plot_nullclines(f,u0,tmax,p;xlims=[-1.0,1.0],ylims=[-1.0,1.0],size=(700,500))
    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    p1 = plot_nullclines(f,p;xlims=xlims,ylims=ylims)
    condition(u,t,integrator) = (u[1]*u[1]+u[2]*u[2]) > max(xrange*xrange,yrange*yrange)
    affect!(integrator) = terminate!(integrator)
    sol = solve(ODEProblem(f,u0,(0.0,tmax),p),callback=DiscreteCallback(condition,affect!))
    plot!(p1,sol,vars=(1,2),c=:black,arrow=true,xlims=xlims,ylims=ylims,xaxis=("μ2"),yaxis=("y"),size=size)
end    

function solve_plot_nullclines_flux(f,tmax,p;Ngrid=10,xlims=[-1.0,1.0],ylims=[-1.0,1.0],size=(700,500))
    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    u0_arr = vec([[xlims[1]+i*xrange/Ngrid,ylims[1]+j*yrange/Ngrid] for i=0:Ngrid, j=0:Ngrid])
    p1 = plot_nullclines(f,p,xlims=xlims,ylims=ylims)
    prob = ODEProblem(f,u0_arr[1],(0.0,tmax),p)
    ensamble_prob = EnsembleProblem(prob,prob_func=(prob,i,repeat;u0=u0_arr)->(remake(prob,u0=u0[i])))
    condition(u,t,integrator) = (u[1]*u[1]+u[2]*u[2]) > max(xrange*xrange,yrange*yrange)
    affect!(integrator) = terminate!(integrator)
    sol = solve(ensamble_prob,EnsembleThreads(),trajectories=length(u0_arr),callback=DiscreteCallback(condition,affect!))
    plot!(p1,sol,vars=(1,2),arrows=true,c=:black,linewidth=0.5,xlims=xlims,ylims=ylims,size=size)
end    

function solve_plot_animated(f,p,N,dt;Ngrid=10,xlims=[-1.0,1.0],ylims=[-1.0,1.0],size=(400,400),nullclines=false)
    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    u0_arr = vec([[xlims[1]+i*xrange/Ngrid,ylims[1]+j*yrange/Ngrid] for i=0:Ngrid, j=0:Ngrid])
    prob = ODEProblem(f,u0_arr[1],(0.0,N*dt),p)
    ensamble_prob = EnsembleProblem(prob,prob_func=(prob,i,repeat;u0=u0_arr)->(remake(prob,u0=u0[i])))
    sol = solve(ensamble_prob,EnsembleThreads(),trajectories=length(u0_arr))
    x = reduce(hcat,[sol[n](0:dt:N*dt,idxs=1).u for n=1:length(sol)])
    y = reduce(hcat,[sol[n](0:dt:N*dt,idxs=2).u for n=1:length(sol)])
    if (nullclines)
        p1 = plot_nullclines(f,p;xlims=xlims,ylims=ylims,size=size) 
    else    
        p1 = plot([x[1:2,:]],[y[1:2,:]],color=:black,xlims=xlims,ylims=ylims,legend=false,size=size)
    end    
    anim = @animate for n=3:N
        plot!(p1,x[n-1:n,:],y[n-1:n,:],color=:black)
    end
    anim
end    

function classification_linear(A;Ngrid=5,tmax=2.0,xlims=[-1.2,1.2],ylims=[-1.2,1.2])
    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    u0_arr = vec([[xlims[1]+i*xrange/Ngrid,ylims[1]+j*yrange/Ngrid] for i=0:Ngrid, j=0:Ngrid])
    prob = ODEProblem((u,p,t) -> A*u, [0;0], (0,tmax))
    ensamble_prob = EnsembleProblem(prob,prob_func=(prob,i,repeat;u0=u0_arr)->(remake(prob,u0=u0[i])))
    sol = solve(ensamble_prob,EnsembleThreads(),trajectories=length(u0_arr))
    p2 = plot(sol,vars=(1,2),xlims=xlims,ylims=ylims,arrows=true,c=:black,linewidth=0.5)
    maxtr = 2.0
    maxdet = 1.3
    trv = -maxtr:maxtr/30:maxtr
    trdet = abs2.(trv)/4
    p1 = plot(trv,trdet,xaxis=("Tr(A)",(-maxtr,maxtr)),yaxis=("Det(A)",(-maxdet,maxdet)),legend=false)
    plot!(p1,[-maxtr,maxtr],[0,0])
    plot!(p1,[0,0],[0,maxdet])
    txt=["Saddle","Stable\nNode","Stable\nFocus","Unstable\nFocus","Unstable\nNode"]
    annotate!([(0,-0.5*maxdet,(txt[1],12,:gray))])
    annotate!([(-0.8*maxtr,0.2*maxdet,(txt[2],10,:gray)),(-0.5*maxtr,0.7*maxdet,(txt[3],10,:gray))])
    annotate!([(0.8*maxtr,0.2*maxdet,(txt[5],10,:gray)),(0.5*maxtr,0.7*maxdet,(txt[4],10,:gray))])
    atxt = string(A[1,1]) * " " * string(A[1,2]) * "\n" * string(A[2,1]) * " " * string(A[2,2])
    scatter!(p1,[tr(A)],[det(A)],c=:red)
    annotate!([(tr(A),det(A)-0.1*maxdet,(atxt,9,:red,:center))])
    plot(p2,p1,layout=(1,2),size=(900,400))
end

function plot_manifolds(p1,f,f_jac,u0_array,p;tmax=30,delta=0.001,repulsor=false,xlims=[-1.0,1.0],ylims=[-1.0,1.0],size=(700,500))
    xrange = xlims[2]-xlims[1]
    yrange = ylims[2]-ylims[1]
    condition(u,t,integrator) = (u[1]*u[1]+u[2]*u[2]) > max(xrange*xrange,yrange*yrange)
    affect!(integrator) = terminate!(integrator)
    for u0 in u0_array
        A = f_jac(u0,p)
        if det(A)<0
            av = eigen(A)
            for n=1:2
                if real(av.values[n])>0
                    u1 = u0+delta*av.vectors[:,n]
                    sol = solve(ODEProblem(f,u1,(0.0,tmax),p),callback=DiscreteCallback(condition,affect!))
                    plot!(p1,sol,vars=(1,2),c=:red)
                    u1 = u0-delta*av.vectors[:,n]
                    sol = solve(ODEProblem(f,u1,(0.0,tmax),p),callback=DiscreteCallback(condition,affect!))
                    plot!(p1,sol,vars=(1,2),c=:red)
                else
                    u1 = u0+delta*av.vectors[:,n]
                    sol = solve(ODEProblem(f,u1,(0.0,-tmax),p),callback=DiscreteCallback(condition,affect!))
                    plot!(p1,sol,vars=(1,2),c=:blue)
                    u1 = u0-delta*av.vectors[:,n]
                    sol = solve(ODEProblem(f,u1,(0.0,-tmax),p),callback=DiscreteCallback(condition,affect!))
                    plot!(p1,sol,vars=(1,2),c=:blue)
                end  
            end
        else
            if tr(A)<0
                scatter!(p1,[u0[1]],[u0[2]],c=:black)
            else
                scatter!(p1,[u0[1]],[u0[2]],c=:red)
                if repulsor
                    u1 = u0.+[delta;delta]
                    sol = solve(ODEProblem(f,u1,(0.0,3*tmax),p),callback=DiscreteCallback(condition,affect!))
                    plot!(p1,sol,vars=(1,2),c=:purple,alpha=0.5)
                end    
            end    
        end    
    end
    plot(p1,legend=false,xlims=xlims,ylims=ylims,size=size)
end

function plot_manifolds(f,f_jac,u0_array,p;tmax=30,delta=0.001,repulsor=false,xlims=[-1.0,1.0],ylims=[-1.0,1.0],size=(700,500))
    p1=plot()
    plot_manifolds(p1,f,f_jac,u0_array,p;tmax=tmax,delta=delta,repulsor=repulsor,xlims=xlims,ylims=ylims,size=size)
end

function phase_portrait(f,p;tmax=50,delta=0.001,xlims=[-1.0,1.0],ylims=[-1.0,1.0],size=(700,500))
    # Find fixedpoints in interval 
    function fsv((x,y))
        du = f(similar([x,y]),[x,y],p,0)
        return SVector(du[1],du[2])
    end    
    X = xlims[1]..xlims[2]
    Y = ylims[1]..ylims[2]
    rts = roots(fsv, X × Y)
    u0_arr=[[mid(rt.interval[1]);mid(rt.interval[2])] for rt in rts]
    p1 = plot_nullclines(f,p;xlims=xlims,ylims=ylims,npts=50,regions=false)
    f_jac(u0,p) = ForwardDiff.jacobian(x -> f(similar(x),x,p,0), u0)
    plot_manifolds(p1,f,f_jac,u0_arr,p;tmax=tmax,delta=delta,repulsor=true,xlims=xlims,ylims=ylims,size=size)
end    
    
