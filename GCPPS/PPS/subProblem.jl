function subProb(a::Vector{Float64},b::Float64,i::Int64,j::Int64,k::Int64,act::Int64,TliS::Number)
    println("\n*************************** Sous problème ************************************************\n")
    pricing=Model(CPLEX.Optimizer)
    xa=@variable(pricing,xa[t in 1:T,p in 1:pow],Bin)
    x=@variable(pricing,x[t in 1:T]>=0)
    E=@variable(pricing,E)
    y=@variable(pricing,y[t in 0:T]>=0) # température
    set_optimizer(pricing, optimizer_with_attributes( CPLEX.Optimizer,"CPX_PARAM_TILIM" =>TliS))
    if act==1
             @constraint(pricing,E<=0.01)
             set_optimizer(pricing, optimizer_with_attributes( CPLEX.Optimizer,"CPX_PARAM_INTSOLLIM"=>1))
    end
    @objective(pricing,Min,E) 
    @constraint(pricing,E==-sum(a[t]*Pa[i,j,k,p]*xa[t,p] for p in 1:pow,t in 1:T)-b )
    @constraint(pricing, [t in 1:T],sum(xa[t,p] for p in 1:pow)<=1)
    @constraint(pricing, [t in 1:T],x[t]==sum(Pa[i,j,k,p]*xa[t,p] for p in 1:pow))
    if j==1
              @constraint(pricing, [t in tf[k,i,u]+1:T],sum(xa[t,p] for p in 1:pow)==0)
              @constraint(pricing, [t in 1:T], y[t]==y[t-1]+(delta/Cm[k,i])*(1000*sum(Pa[i,1,k,p]*xa[t,p] for p in 1:pow)-U[k,i]*(y[t-1]-Yout[t]))) # contrainte température maison
              @constraint(pricing, y[0]==y0[k,i])
              for l in 1:u
                        @constraint(pricing, [t in td[k,i,l]:tf[k,i,l]], tmin[k,i,l]<=y[t]<= tmax[k,i,l])
             end
   else
             @constraint(pricing, y[0]==ti[k,i])      
             @constraint(pricing,[l in 1:u], y[tt+48*(l-1)]==ti[k,i])
             @constraint(pricing, [l in 1:u],temp_min<=y[48*(l-1)+ec]<= temp_max)
             @constraint(pricing,[t in 1:T;t ∉ [tt+48*(l-1) for l in 1:u]], y[t]==(y[t-1]*r*cpo*Masse[k,i]+S[k,i]*K*dist*r*Ta+ sum(Pa[i,2,k,p]*xa[t,p] for p in 1:pow)*dist*860)/(Masse[k,i]*cpo*r+S[k,i]*dist*K*r)) # température eau 

   end      
   optimize!(pricing)

  return JuMP.value.(x),JuMP.value.(E), JuMP.objective_value(pricing),round(MOI.get(pricing, MOI.SolveTimeSec()),digits=2)
end
  
