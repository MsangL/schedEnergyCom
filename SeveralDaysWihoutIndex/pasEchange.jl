function PasEchange(jour::Int,Tli::Number)
    @info "\n********************************* When members do not exchange ***********************************************\n"
    mod=Model(CPLEX.Optimizer)
    set_optimizer(mod, optimizer_with_attributes( CPLEX.Optimizer,"CPX_PARAM_TILIM" =>Tli))
    Gain=@variable(mod,Gain[i in 1:n])
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
    inj=@variable(mod,inj[t  in 1:T]>=0)
    ext=@variable(mod,ext[t  in 1:T]>=0)
    Puis=@variable(mod,Puis[i in 1:n,j in 1:JA,k in 1:Nbp, t in 1:T]>=0)
    tot=@variable(mod,tot[t in 1:T]>=0)
    Bess=@variable(mod,Bess[i in 1:T]>=0)
  ############
    @objective(mod,Min,dist*sum(Cg[i,t]+Pg[i,t] for i in 1:n,t in 1:T))
    @constraint(mod,[i in 1:n, t in 1:T],sum(Pa[i,j,k,s]*xa[j,i,k,t,s] for j in 1:JA,s in 1:pow, k in 1:Nbp)+sum(Pb[i,t,j,s]*x[i,j,s] for j in 1:JB,s in 1:M )+sum((q[b,i,t]/dist) for b in 1:B)+p[i,t]==L*prod1[i,t]+sum(f[j,i,t] for j in 1:n if i!=j)+Cg[i,t]-Pg[i,t]-sum(f[i,j,t] for j in 1:n if i!=j))
    @constraint(mod,[i in 1:n,t in 1:T],sum(f[i,j,t] for j in 1:n if j!=i)+Pg[i,t]<=(L*prod1[i,t]+sum(E[b,i,t]/dist for b in 1:B)*omega[i]))
    @constraint(mod,[i in 1:n,j in 1:JB],sum(x[i,j,s] for s in 1:M)==1) #b
    @constraint(mod,[ j in 1:JA, i in 1:n , t in 1:T,k in 1:Nbp], sum(xa[j,i,k,t,s] for s in 1:pow)<=1) #
    
    @constraint(mod, [ t in 1:T, i in 1:n,k in 1:Nbp], y[i,k,t,1]==y[i,k,t-1,1]+(delta/Cm[k,i])*(1000*sum(Pa[i,1,k,s]*xa[1,i,k,t,s]*JobA[i,1,k] for s in 1:pow)-U[k,i]*(y[i,k,t-1,1]-Yout[t]))) # contrainte température maison
   @constraint(mod, [ i in 1:n, k in 1:Nbp,t in 1:T; t!=tt], y[i,k,t,2]==(y[i,k,t-1,2]*r*Masse[k,i]+S[k,i]*K*dist*r*Ta+ sum(Pa[i,2,k,s]*xa[2,i,k,t,s]*JobA[i,2,k] for s in 1:pow)*dist*860)/(Masse[k,i]*cpo*r+S[k,i]*dist*K*r)) # température eau 
  
   @constraint(mod,[ i in 1:n, k in 1:Nbp], y[i,k,tt,2]==ti[k,i])
   @constraint(mod,[ i in 1:n,k in 1:Nbp], y[i,k,0,1]==y0[k,i])
    if jour==1
                  @constraint(mod,[ i in 1:n, k in 1:Nbp], y[i,k,0,2]==ti[k,i])
     else  
                  @constraint(mod,[ i in 1:n, k in 1:Nbp], y[i,k,0,2]==te0[k,i])

     end
     
    @constraint(mod,[ i in 1:n,k in 1:Nbp ; JobA[i,2,k]>0], temp_min<=y[i,k,ec,2]<= temp_max)
    @constraint(mod,[ i in 1:n,k in 1:Nbp, t in td[k,i]:tf[k,i]; JobA[i,1,k]>0], tmin[k,i]<=y[i,k,t,1]<= tmax[k,i])
    @constraint(mod,[i in 1:n,j in 1:JA,k in 1:Nbp],sum(xa[j,i,k,t,p] for t in 1:T,p in 1:pow)<=T*JobA[i,j,k])           # en parler avc Michael
    #####
    @constraint(mod,[i in 1:n,b in 1:B,t in 1:T],dist*d*Pd[i,b]*(z[b,i,t]-1)<=q[b,i,t])
    @constraint(mod,[i in 1:n,b in 1:B,t in 1:T],q[b,i,t]<=dist*c*Pc[i,b]*z[b,i,t])
    @constraint(mod,[i in 1:n,b in 1:B,t in 1:T-1],z[b,i,t+1]-z[b,i,t]<=w[b,i,t])
    @constraint(mod,[i in 1:n,b in 1:B],sum(w[b,i,t] for t in 1:T)<=fi)
    
    @constraint(mod,[i in 1:n,b in 1:B,t in 1:T],E[b,i,t]==eta*E[b,i,t-1]+q[b,i,t])
    @constraint(mod,[i in 1:n,b in 1:B],E[b,i,0]==xi0[i,b])
    @constraint(mod,[i in 1:n,t in 1:T,b in 1:B],E[b,i,t]<=gamma[i,b])
    
    @constraint(mod,[i in 1:n,t in 1:T],sum(f[i,j,t] for j in 1:n if i!=j )+Pg[i,t]<=P_sous[i])
    @constraint(mod,[i in 1:n,t in 1:T],sum(f[j,i,t] for j in 1:n if i!=j )+Cg[i,t]<=P_sous[i])
    @constraint(mod,[j in 1:n,i in 1:n,t in 1:T],f[j,i,t]==0)                                                # f_ij=0
    @constraint(mod,[i in 1:n],Gain[i]==(Pv*V[i]-Pedf*C[i]))                                           # gain avant
    ##########################################################################################ontraintes auxilliaires
    @constraint(mod,[i in 1:n],V[i]==dist*sum(Pg[i,t] for t in 1:T ))
    @constraint(mod,[i in 1:n],C[i]==dist*sum(Cg[i,t] for t in 1:T ))
    @constraint(mod,[j in 1:JA,t in 1:T, i in 1:n,k in 1:Nbp],Puis[i,j,k,t]==sum(Pa[i,j,k,s]*xa[j,i,k,t,s] for s in 1:pow))
    @constraint(mod,[t in 1:T],tot[t]==dist*sum(Puis[i,j,k,t] for i in 1:n, j in 1:JA,k in 1:Nbp)+dist*sum(p[i,t] for i in 1:n)+dist*sum(Pb[i,t,j,s] for i in 1:n,j in 1:JB,s in 1:M if JB>1))
    @constraint(mod,[t in 1:T],inj[t]==dist*sum(Pg[i,t] for i in 1:n ))
    @constraint(mod,[t in 1:T],ext[t]==dist*sum(Cg[i,t] for i in 1:n ))
    @constraint(mod,[t in 1:T],Bess[t]==sum(E[b,i,t] for i in 1:n, b in B))
    optimize!(mod)  
    return (JuMP.value.(C),JuMP.value.(Gain),JuMP.value.(E),JuMP.value.(y),MOI.get(mod, MOI.RelativeGap()),JuMP.value.(inj),JuMP.value.(ext),JuMP.value.(Bess),JuMP.value.(tot),JuMP.value.(q))
end
