function MasterProb(P,Ma,GainA,state,R,Tli)
    mod=Model(CPLEX.Optimizer)
    set_optimizer(mod, optimizer_with_attributes( CPLEX.Optimizer,"CPX_PARAM_TILIM" =>Tli)) #3600
    set_silent(mod)

    Gain=@variable(mod,Gain[i in 1:n])
     if state==1
                   println("\n***************************************** Master problrem with integrality constraints ***************************************\n")
                   # this one is solved after the last iteration of the column generation
                   xa=@variable(mod,xa[i in 1:n,j in 1:JA,k in 1:Nbp,s in 1:Ma]>=0)
                   z=@variable(mod,z[b in 1:B,i in 1:n,t in 1:T]>=0,Bin)
                   x=@variable(mod,x[i in 1:n,j in 1:JB,s in 1:M]>=0,Bin)
                   w=@variable(mod,w[b in 1:B,i in 1:n,t in 1:T]>=0,Bin)
                   @constraint(mod,bet[i in 1:n,j in 1:JA,k in 1:Nbp; sum(JobA[i,j,k,l] for l in 1:u)>=1],sum(xa[i,j,k,s] for s in 1:Ma) >= 1)
    else
                   println("\n***************************************** Master problem ***************************************\n")
                   xa=@variable(mod,xa[i in 1:n,j in 1:JA,k in 1:Nbp,s in 1:Ma])
                   z=@variable(mod,z[b in 1:B,i in 1:n,t in 1:T])
                   x=@variable(mod,x[i in 1:n,j in 1:JB,s in 1:M])
                   w=@variable(mod,w[b in 1:B,i in 1:n,t in 1:T])
                   @constraint(mod,[i in 1:n, j in 1:JA,k in 1:Nbp,s in 1:Ma],0<=xa[i,j,k,s]<=1)
                   @constraint(mod,[i in 1:n, j in 1:JB,s in 1:M],0<=x[i,j,s]<=1)
                   @constraint(mod,[i in 1:n, b in 1:B,t in 1:T],0<=z[b,i,t]<=1)
                   @constraint(mod,[i in 1:n, b in 1:B,t in 1:T],0<=w[b,i,t]<=1)
                   @constraint(mod,bet[i in 1:n,j in 1:JA,k in 1:Nbp; sum(JobA[i,j,k,l] for l in 1:u)>=1],sum(xa[i,j,k,s] for s in 1:Ma) == 1)

    end

    q=@variable(mod,q[b in 1:B,i in 1:n,t in 1:T])
    f=@variable(mod,f[i in 1:n,j in 1:n,t in 1:T]>=0)
    Pg=@variable(mod,Pg[i in 1:n, t  in 1:T]>=0)
    Cg=@variable(mod,Cg[i in 1:n, t  in 1:T]>=0)
    E=@variable(mod,E[b in 1:B,i in 1:n,t in 0:T]>=0)
    Fsend=@variable(mod,Fsend[i in 1:n]>=0)
    Freceive=@variable(mod,Freceive[i in 1:n]>=0)
    C=@variable(mod,C[i in 1:n]>=0)
    V=@variable(mod,V[i in 1:n]>=0)
    ####################################
    inj=@variable(mod,inj[t  in 1:T]>=0)
    ext=@variable(mod,ext[t  in 1:T]>=0)
    Puis=@variable(mod,Puis[i in 1:n,j in 1:JA,k in 1:Nbp, t in 1:T]>=0)
    tot=@variable(mod,tot[t in 1:T]>=0)
    Bess=@variable(mod,Bess[i in 1:T]>=0)
    ####################################
    @objective(mod,Min,dist*sum(Cg[i,t] for i in 1:n,t in 1:T))
    @constraint(mod,alpha[i in 1:n, t in 1:T],sum(xa[i,j,k,s]*P[i,j,k,t,s] for j in 1:JA,s in 1:Ma, k in 1:Nbp)+sum(Pb[i,t,j,s]*x[i,j,s] for j in 1:JB,s in 1:M if JB>0)+sum((q[b,i,t]/dist) for b in 1:B)+p[i,t]==L*prod1[i,t]+sum(f[j,i,t] for j in 1:n if i!=j)+Cg[i,t]-Pg[i,t]-sum(f[i,j,t] for j in 1:n if i!=j))

    @constraint(mod,[i in 1:n,j in 1:JB],sum(x[i,j,s] for s in 1:M)>=1) 
    @constraint(mod,[i in 1:n,t in 1:T],sum(f[i,j,t] for j in 1:n if j!=i)+Pg[i,t]<=(L*prod1[i,t]+sum(E[b,i,t]/dist for b in 1:B)*omega[i]))
    @constraint(mod,[i in 1:n,b in 1:B,t in 1:T],dist*d*Pd[i,b]*(z[b,i,t]-1)<=q[b,i,t])
    @constraint(mod,[i in 1:n,b in 1:B,t in 1:T],q[b,i,t]<=dist*c*Pc[i,b]*z[b,i,t])
    @constraint(mod,[i in 1:n,b in 1:B,t in 1:T-1],z[b,i,t+1]-z[b,i,t]<=w[b,i,t])
    @constraint(mod,[i in 1:n,b in 1:B],sum(w[b,i,t] for t in 1:T)<=fi)
    @constraint(mod,[i in 1:n,b in 1:B,t in 1:T],E[b,i,t]==eta*E[b,i,t-1]+q[b,i,t])
    @constraint(mod,[i in 1:n,b in 1:B],E[b,i,0]==xi[i,b])
    @constraint(mod,[i in 1:n,b in 1:B],E[b,i,T]==xi[i,b])
    @constraint(mod,[i in 1:n,t in 1:T,b in 1:B],E[b,i,t]<=gamma[i,b])
    @constraint(mod,[i in 1:n],Fsend[i]==dist*sum(f[i,j,t] for j in 1:n,t in 1:T if i!=j))
    @constraint(mod,[i in 1:n],Freceive[i]==dist*sum(f[j,i,t] for j in 1:n,t in 1:T if i!=j))
    @constraint(mod,[i in 1:n],V[i]==dist*sum(Pg[i,t] for t in 1:T ))
    @constraint(mod,[i in 1:n],C[i]==dist*sum(Cg[i,t] for t in 1:T ))
    @constraint(mod,[i in 1:n,t in 1:T],sum(f[i,j,t] for j in 1:n if i!=j )+Pg[i,t]<=P_sous[i])
    @constraint(mod,[i in 1:n,t in 1:T],sum(f[j,i,t] for j in 1:n if i!=j )+Cg[i,t]<=P_sous[i])
    @constraint(mod,[i in 1:n],Gain[i]==(Pvc*Fsend[i]-Pac*Freceive[i]+oui*V[i]-Pedf*C[i]))    # Gain après
                           @constraint(mod,[i in 1:n],Gain[i]>=GainA[i]-beta*abs(GainA[i]))   #gain après sup a 15% gain avant

    ################################# 
    @constraint(mod,[j in 1:JA,t in 1:T, i in 1:n,k in 1:Nbp],Puis[i,j,k,t]==sum(P[i,j,k,t,s]*xa[i,j,k,s] for s in 1:Ma)) 
    @constraint(mod,[t in 1:T],tot[t]==dist*sum(xa[i,j,k,s]*P[i,j,k,t,s] for i in 1:n, j in 1:JA,s in 1:Ma, k in 1:Nbp if Nbp>0)+dist*sum(p[i,t] for i in 1:n)+dist*sum(Pb[i,t,j,s]*x[i,j,s] for i in 1:n,j in 1:JB,s in 1:M if JB>1))
    @constraint(mod,[t in 1:T],inj[t]==dist*sum(Pg[i,t] for i in 1:n ))
    @constraint(mod,[t in 1:T],ext[t]==dist*sum(Cg[i,t] for i in 1:n ))
    @constraint(mod,[t in 1:T],Bess[t]==sum(E[b,i,t] for i in 1:n, b in B))
    optimize!(mod)    
  
    if state!=1
                    return (mod,alpha,bet,JuMP.objective_value(mod),round(MOI.get(mod, MOI.SolveTimeSec()),digits=2))
    else 
                   return (JuMP.objective_value(mod),JuMP.value.(C),JuMP.value.(V),JuMP.value.(Gain),JuMP.value.(tot),JuMP.value.(inj),JuMP.value.(ext),JuMP.value.(Bess),JuMP.value.(x),JuMP.value.(Puis), round(MOI.get(mod, MOI.SolveTimeSec()),digits=2), JuMP.value.(q))
    end
end


    
 
