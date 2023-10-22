#function subProb(a,b,act,R,TliS)
function subProb(a,b,act::Number,R::Array{Float64},TliS::Number) 
    # This function returns the first feasible solution if act equals to 1
    pricing=Model(CPLEX.Optimizer)
    x=@variable(pricing,x[i in 1:n,j in 1:JA,k in 1:Nbp,t in 1:T]>=0)
    E=@variable(pricing,E[i in 1:n,j in 1:JA,k in 1:Nbp])
    y=@variable(pricing,y[i in 1:n,k in 1:Nbp, t in 0:T, j in 1:JA]>=0) # température
    set_optimizer(pricing, optimizer_with_attributes( CPLEX.Optimizer,"CPX_PARAM_TILIM" =>TliS))
    psi=@variable(pricing,psi[i in 1:n,j in 1:JA,k in 1:Nbp,t in 1:T,p in 1:pow], Bin)


      if act==1
             @constraint(pricing,[i in 1:n,j in 1:JA,k in 1:Nbp; R[i,j,k]>0],E[i,j,k]<=0.01)
             set_optimizer(pricing, optimizer_with_attributes( CPLEX.Optimizer,"CPX_PARAM_INTSOLLIM"=>1))
      end

                @objective(pricing,Min, sum(E[i,j,k] for i in 1:n,j in 1:JA,k in 1:Nbp if R[i,j,k]>0)) 

              @constraint(pricing,[ i in 1:n,j in 1:JA,k in 1:Nbp;R[i,j,k]>0],E[i,j,k]==-sum(dual(a[i,t])*Pa[i,j,k,p]*(psi[i,j,k,t,p]-psi[i,j,k,t,p+1]) for p in 1:pow-1,t in 1:T)-sum(dual(a[i,t])*Pa[i,j,k,pow]*psi[i,j,k,t,pow] for t in 1:T)-dual(b[i,j,k]) )
              @constraint(pricing, [i in 1:n,j in 1:JA,t in 1:T,k in 1:Nbp],x[i,j,k,t]==sum(Pa[i,j,k,p]*(psi[i,j,k,t,p]-psi[i,j,k,t,p+1]) for p in 1:pow-1)+Pa[i,j,k,pow]*psi[i,j,k,t,pow])

              #Fortz Mars 16th, 2022
              @constraint(pricing,[i in 1:n,j in 1:JA,k in 1:Nbp, t in 1:T, p in 2:pow],psi[i,j,k,t,p]<=psi[i,j,k,t,p-1])
              @constraint(pricing,[i in 1:n,j in 1:JA,k in 1:Nbp, t in 1:T; R[i,j,k]>0],psi[i,j,k,t,1]==1)
               
              @constraint(pricing, [t in 1:T,i in 1:n,k in 1:Nbp], y[i,k,t,1]==y[i,k,t-1,1]+(delta/Cm[k,i])*(1000*(sum(Pa[i,1,k,p]*(psi[i,1,k,t,p]-psi[i,1,k,t,p+1]) for p in 1:pow-1)+Pa[i,1,k,pow]*psi[i,1,k,t,pow])-U[k,i]*(y[i,k,t-1,1]-Yout[t]))) 
              @constraint(pricing, [i in 1:n,k in 1:Nbp], y[i,k,0,1]==y0[k,i])
                 
              @constraint(pricing, [i in 1:n,k in 1:Nbp], y[i,k,0,2]==ti[k,i])
              @constraint(pricing, [i in 1:n,k in 1:Nbp,l in 1:u], y[i,k,tt+48*(l-1),2]==ti[k,i])
              @constraint(pricing, [i in 1:n,k in 1:Nbp, t in 1:T;   t ∉ [tt+48*(l-1) for l in 1:u]  ], y[i,k,t,2]==(y[i,k,t-1,2]*r*cpo*Masse[k,i]+S[k,i]*K*dist*r*Ta+(sum(Pa[i,2,k,p]*(psi[i,2,k,t,p]-psi[i,2,k,t,p+1]) for p in 1:pow-1)+Pa[i,2,k,pow]*psi[i,2,k,t,pow])*dist*860)/(Masse[k,i]*cpo*r+S[k,i]*dist*K*r)) 
                
              @constraint(pricing, [i in 1:n,k in 1:Nbp,l in 1:u; JobA[i,2,k,l]>0 ], temp_min<=y[i,k,48*(l-1)+ec,2]<= temp_max)
              @constraint(pricing, [i in 1:n,k in 1:Nbp,l in 1:u ,t in td[k,i,l]:tf[k,i,l]; JobA[i,1,k,l]>0], tmin[k,i,l]<=y[i,k,t,1]<= tmax[k,i,l])

    		optimize!(pricing)
              return JuMP.value.(x),JuMP.value.(E), JuMP.objective_value(pricing),round(MOI.get(pricing, MOI.SolveTimeSec()),digits=2),MOI.get(pricing, MOI.RelativeGap())
  end
