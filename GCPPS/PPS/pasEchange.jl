function pasEchange(Tlimit::Number)
      mod=Model(CPLEX.Optimizer)
      set_optimizer(mod, optimizer_with_attributes( CPLEX.Optimizer,"CPX_PARAM_TILIM" =>Tlimit))
      Gain=@variable(mod,Gain[i in 1:n])
      #set_silent(mod)
      x=@variable(mod,x[i in 1:n,j in 1:JB,s in 1:M]>=0,Bin)
      xa=@variable(mod,xa[j in 1:JA, i in 1:n,k in 1:Nbp,t in 1:T,s in 1:pow],Bin)
      z=@variable(mod,z[b in 1:B,i in 1:n,t in 1:T]>=0,Bin)
      q=@variable(mod,q[b in 1:B,i in 1:n,t in 1:T])
      w=@variable(mod,w[b in 1:B,i in 1:n,t in 1:T]>=0,Bin)
      f=@variable(mod,f[i in 1:n,j in 1:n,t in 1:T]>=0)
      Pg=@variable(mod,Pg[i in 1:n, t  in 1:T]>=0)
      Cg=@variable(mod,Cg[i in 1:n, t  in 1:T]>=0)
      E=@variable(mod,E[b in 1:B,i in 1:n,t in 0:T]>=0)
      y=@variable(mod,y[i in 1:n,k in 1:Nbp, t in 0:T, j in 1:JA]>=0) # température 
      C=@variable(mod,C[i in 1:n]>=0)
      V=@variable(mod,V[i in 1:n]>=0)
     
      ############
      @objective(mod,Min,dist*sum(Cg[i,t] for i in 1:n,t in 1:T))
      @constraint(mod,[i in 1:n, t in 1:T],sum(Pa[i,j,k,s]*xa[j,i,k,t,s] for j in 1:JA,s in 1:pow, k in 1:Nbp)+sum(Pb[i,t,j,s]*x[i,j,s] for j in 1:JB,s in 1:M )+sum((q[b,i,t]/dist) for b in 1:B)+p[i,t]==L*prod1[i,t]+sum(f[j,i,t] for j in 1:n if i!=j)+Cg[i,t]-Pg[i,t]-sum(f[i,j,t] for j in 1:n if i!=j))
      @constraint(mod,[i in 1:n,t in 1:T],sum(f[i,j,t] for j in 1:n if j!=i)+Pg[i,t]<=(L*prod1[i,t]+sum(E[b,i,t]/dist for b in 1:B)*omega[i]))
      @constraint(mod,[i in 1:n,j in 1:JB],sum(x[i,j,s] for s in 1:M)==1) #b
      @constraint(mod,[ j in 1:JA, i in 1:n , t in 1:T,k in 1:Nbp], sum(xa[j,i,k,t,s] for s in 1:pow)<=1) #
      
      @constraint(mod, [ t in 1:T, i in 1:n,k in 1:Nbp ], y[i,k,t,1]==y[i,k,t-1,1]+(delta/Cm[k,i])*(1000*sum(Pa[i,1,k,s]*xa[1,i,k,t,s] for s in 1:pow)-U[k,i]*(y[i,k,t-1,1]-Yout[t]))) # contrainte température maison
      @constraint(mod,[ i in 1:n,k in 1:Nbp], y[i,k,0,1]==y0[k,i])
      @constraint(mod,[ i in 1:n, k in 1:Nbp], y[i,k,0,2]==ti[k,i])

      @constraint(mod,[ i in 1:n,k in 1:Nbp,l in 1:u; JobA[i,2,k,l]>0], temp_min<=y[i,k,ec+(l-1)*48,2]<= temp_max)
      @constraint(mod,[ i in 1:n,k in 1:Nbp,l in 1:u, t in td[k,i,l]:tf[k,i,l]; JobA[i,1,k,l]>0], tmin[k,i,l]<=y[i,k,t,1]<= tmax[k,i,l])

     @constraint(mod,[ i in 1:n, k in 1:Nbp,l in 1:u], y[i,k,tt+48*(l-1),2]==ti[k,i])
     @constraint(mod,[ i in 1:n,k in 1:Nbp, t in 1:T; t ∉ [tt+48*(l-1) for l in 1:u] ], y[i,k,t,2]==(y[i,k,t-1,2]*r*cpo*Masse[k,i]+S[k,i]*K*dist*r*Ta+sum(Pa[i,2,k,s]*xa[2,i,k,t,s] for s in 1:pow)*dist*860)/(Masse[k,i]*cpo*r+S[k,i]*dist*K*r)) 
                            
      @constraint(mod,[i in 1:n,j in 1:JA,k in 1:Nbp],sum(xa[j,i,k,t,p] for t in 1:T,p in 1:pow)<=T*sum(JobA[i,j,k,l] for l in 1:u))
      @constraint(mod,[i in 1:n,b in 1:B,t in 1:T],dist*d*Pd[i,b]*(z[b,i,t]-1)<=q[b,i,t])
      @constraint(mod,[i in 1:n,b in 1:B,t in 1:T],q[b,i,t]<=dist*c*Pc[i,b]*z[b,i,t])
      @constraint(mod,[i in 1:n,b in 1:B,t in 1:T-1],z[b,i,t+1]-z[b,i,t]<=w[b,i,t])
      @constraint(mod,[i in 1:n,b in 1:B],sum(w[b,i,t] for t in 1:T)<=fi)
      @constraint(mod,[i in 1:n,b in 1:B,t in 1:T],E[b,i,t]==eta*E[b,i,t-1]+q[b,i,t])
      @constraint(mod,[i in 1:n,b in 1:B],E[b,i,0]==xi[i,b])
      @constraint(mod,[i in 1:n,b in 1:B],E[b,i,T]==xi[i,b])
      @constraint(mod,[i in 1:n,t in 1:T,b in 1:B],E[b,i,t]<=gamma[i,b])
      @constraint(mod,[i in 1:n,t in 1:T],sum(f[i,j,t] for j in 1:n if i!=j )+Pg[i,t]<=P_sous[i])
      @constraint(mod,[i in 1:n,t in 1:T],sum(f[j,i,t] for j in 1:n if i!=j )+Cg[i,t]<=P_sous[i])
      @constraint(mod,[j in 1:n,i in 1:n,t in 1:T],f[j,i,t]==0)                                                # f_ij=0
      @constraint(mod,[i in 1:n],Gain[i]==(Pv*V[i]-Pedf*C[i]))                                           # gain avant
      ##########################################################################################contraintes auxilliaires
      @constraint(mod,[i in 1:n],V[i]==dist*sum(Pg[i,t] for t in 1:T ))
      @constraint(mod,[i in 1:n],C[i]==dist*sum(Cg[i,t] for t in 1:T ))
      
      optimize!(mod)
      return JuMP.value.(Gain),JuMP.value.(xa)
end

