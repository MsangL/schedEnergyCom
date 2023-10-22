function subProb(a,b,act::Number,R::Array{Float64},TliS::Number) 
# in this function, if act equals one, we return the first feasible solution found. It is equal to one if there no such condition
    println("\n*************************** Subproblem ************************************************\n")
    pricing=Model(CPLEX.Optimizer)
    #set_silent(pricing)
    xa=@variable(pricing,xa[j in 1:JA, i in 1:n,k in 1:Nbp,t in 1:T,p in 1:pow],Bin)
    x=@variable(pricing,x[i in 1:n,j in 1:JA,k in 1:Nbp,t in 1:T]>=0)
    E=@variable(pricing,E[i in 1:n,j in 1:JA,k in 1:Nbp])
    y=@variable(pricing,y[i in 1:n,k in 1:Nbp, t in 0:T, j in 1:JA]>=0) # température
    set_optimizer(pricing, optimizer_with_attributes( CPLEX.Optimizer,"CPX_PARAM_TILIM" =>TliS))
    
      if act==1
        @constraint(pricing,[i in 1:n,j in 1:JA,k in 1:Nbp; R[i,j,k]>0],E[i,j,k]<=0.01)
        set_optimizer(pricing, optimizer_with_attributes( CPLEX.Optimizer,"CPX_PARAM_INTSOLLIM"=>1))
      end
    
      @objective(pricing,Min,sum(E[i,j,k] for i in 1:n,j in 1:JA,k in 1:Nbp if R[i,j,k]>0)) 
      @constraint(pricing,[ i in 1:n,j in 1:JA,k in 1:Nbp;R[i,j,k]>0],E[i,j,k]==-sum(dual(a[i,t])*Pa[i,j,k,p]*xa[j,i,k,t,p] for p in 1:pow,t in 1:T)-dual(b[i,j,k]) )
    
      @constraint(pricing, [i in 1:n,j in 1:JA,t in 1:T,k in 1:Nbp;R[i,j,k]>0],sum(xa[j,i,k,t,p] for p in 1:pow)<=1)
      @constraint(pricing, [i in 1:n,j in 1:JA,t in 1:T,k in 1:Nbp],x[i,j,k,t]==sum(Pa[i,j,k,p]*xa[j,i,k,t,p] for p in 1:pow))
    
      @constraint(pricing, [t in 1:T,i in 1:n,k in 1:Nbp], y[i,k,t,1]==y[i,k,t-1,1]+(delta/Cm[k,i])*(1000*sum(Pa[i,1,k,p]*xa[1,i,k,t,p] for p in 1:pow)-U[k,i]*(y[i,k,t-1,1]-Yout[t]))) # contrainte température maison
      @constraint(pricing, [i in 1:n,k in 1:Nbp], y[i,k,0,1]==y0[k,i])
      @constraint(pricing, [i in 1:n,k in 1:Nbp], y[i,k,0,2]==ti[k,i])
      @constraint(pricing, [i in 1:n,k in 1:Nbp,l in 1:u], y[i,k,tt+48*(l-1),2]==ti[k,i])
      @constraint(pricing, [i in 1:n,k in 1:Nbp, t in 1:T;   t ∉ [tt+48*(l-1) for l in 1:u]  ], y[i,k,t,2]==(y[i,k,t-1,2]*r*cpo*Masse[k,i]+S[k,i]*K*dist*r*Ta+ sum(Pa[i,2,k,p]*xa[2,i,k,t,p] for p in 1:pow)*dist*860)/(Masse[k,i]*cpo*r+S[k,i]*dist*K*r)) 
                
      @constraint(pricing, [i in 1:n,k in 1:Nbp,l in 1:u; JobA[i,2,k,l]>0 ], temp_min<=y[i,k,48*(l-1)+ec,2]<= temp_max)
      @constraint(pricing, [i in 1:n,k in 1:Nbp,l in 1:u ,t in td[k,i,l]:tf[k,i,l]; JobA[i,1,k,l]>0], tmin[k,i,l]<=y[i,k,t,1]<= tmax[k,i,l])

    optimize!(pricing)
    return JuMP.value.(x),JuMP.value.(E), JuMP.objective_value(pricing),round(MOI.get(pricing, MOI.SolveTimeSec()),digits=2),MOI.get(pricing, MOI.RelativeGap())
    
    end
