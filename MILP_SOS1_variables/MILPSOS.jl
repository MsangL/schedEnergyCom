### By SANGARE Mariam ###
### November 2021

 using JuMP,CPLEX
 
include("pasEchange.jl")

function AvecEchange(GainA::Vector{Float64},jour::Int,Tli::Number)
    println("\n***************************************** Members exchange their surplus ***************************************\n")
      mod=Model(CPLEX.Optimizer)
      set_optimizer(mod, optimizer_with_attributes( CPLEX.Optimizer,"CPX_PARAM_TILIM" =>Tli, "CPX_PARAM_WORKMEM"=>500,"CPX_PARAM_WORKDIR"=> "Set this param","CPX_PARAM_NODEFILEIND"=>3))
      Gain=@variable(mod,Gain[i in 1:n,l in 1:u])
      x=@variable(mod,x[i in 1:n,j in 1:JB,s in 1:M,l in 1:u]>=0,Bin)
      z=@variable(mod,z[b in 1:B,i in 1:n,t in 1:T,l in 1:u]>=0,Bin)
      q=@variable(mod,q[b in 1:B,i in 1:n,t in 1:T,l in 1:u])
      w=@variable(mod,w[b in 1:B,i in 1:n,t in 1:T,l in 1:u]>=0,Bin)
      f=@variable(mod,f[i in 1:n,j in 1:n,t in 1:T,l in 1:u]>=0)
      Pg=@variable(mod,Pg[i in 1:n, t  in 1:T,l in 1:u]>=0)
      Cg=@variable(mod,Cg[i in 1:n, t  in 1:T,l in 1:u]>=0)
      E=@variable(mod,E[b in 1:B,i in 1:n,t in 0:T,l in 1:u]>=0)
      y=@variable(mod,y[i in 1:n,k in 1:Nbp, t in 0:T, j in 1:JA,l in 1:u]>=0) # température 
      C=@variable(mod,C[i in 1:n, l in 1:u]>=0)
      V=@variable(mod,V[i in 1:n, l in 1:u]>=0)
      psi=@variable(mod,psi[i in 1:n,j in 1:JA,k in 1:Nbp,t in 1:T,p in 1:pow, l in 1:u], Bin)
      Fsend=@variable(mod,Fsend[i in 1:n,l in 1:u]>=0)
      Freceive=@variable(mod,Freceive[i in 1:n,l in 1:u]>=0)

      @objective(mod,Min,dist*sum(Cg[i,t,l] for i in 1:n,t in 1:T, l in 1:u))
      @constraint(mod,[i in 1:n, t in 1:T,l in 1:u],sum(Pa[i,j,k,s]*(psi[i,j,k,t,s,l]-psi[i,j,k,t,s+1,l]) for j in 1:JA,s in 1:pow-1, k in 1:Nbp if JobA[i,j,k,l]>0)+sum(Pa[i,j,k,pow]*psi[i,j,k,t,pow,l] for j in 1:JA, k in 1:Nbp)+sum(Pb[i,t,j,s,l]*x[i,j,s,l] for j in 1:JB,s in 1:M )+sum((q[b,i,t,l]/dist) for b in 1:B)+p[i,t,l]==L*prod1[i,t,l]+sum(f[j,i,t,l] for j in 1:n if i!=j)+Cg[i,t,l]-Pg[i,t,l]-sum(f[i,j,t,l] for j in 1:n if i!=j))
      
      @constraint(mod,[i in 1:n,t in 1:T,l in 1:u],sum(f[i,j,t,l] for j in 1:n if j!=i)+Pg[i,t,l]<=(L*prod1[i,t,l]+sum(E[b,i,t,l]/dist for b in 1:B)*omega[i]))
      @constraint(mod,[i in 1:n,j in 1:JB,l in 1:u],sum(x[i,j,s,l] for s in 1:M)==1) #b

      @constraint(mod,[ i in 1:n,k in 1:Nbp; JobA[i,1,k,1]>0], y[i,k,0,1,1]==y0[k,i])
     @constraint(mod,[ i in 1:n, k in 1:Nbp;JobA[i,2,k,1]>0], y[i,k,0,2,1]==ti[k,i,1])
      @constraint(mod,[ i in 1:n,k in 1:Nbp,l in 2:u; u>1], y[i,k,0,1,l]==y[i,k,T,1,l-1])
      @constraint(mod,[ i in 1:n, k in 1:Nbp,l in 2:u; u>1], y[i,k,0,2,l]==y[i,k,T,2,l-1])
      
      #Fortz Mars 16th, 2022
      @constraint(mod,[i in 1:n,j in 1:JA,k in 1:Nbp, t in 1:T, p in 2:pow, l in 1:u],psi[i,j,k,t,p,l]<=psi[i,j,k,t,p-1,l])
      @constraint(mod,[i in 1:n,j in 1:JA,k in 1:Nbp, t in 1:T, l in 1:u; JobA[i,j,k,l]>0],psi[i,j,k,t,1,l]==1)
      
      @constraint(mod, [ t in 1:T, i in 1:n,k in 1:Nbp,l in 1:u], y[i,k,t,1,l]==y[i,k,t-1,1,l]+(delta/Cm[k,i])*(1000*(sum((Pa[i,1,k,p]-Pa[i,1,k,p-1])*psi[i,1,k,t,p,l] for p in 2:pow )+Pa[i,1,k,1]*psi[i,1,k,t,1,l])-U[k,i]*(y[i,k,t-1,1,l]-Yout[l,t]))) 
      @constraint(mod, [ i in 1:n, k in 1:Nbp,l in 1:u,  t in 1:T;t!=tt], y[i,k,t,2,l]==(y[i,k,t-1,2,l]*r*cpo*Masse[k,i]+S[k,i]*K*dist*r*Ta+ ( sum((Pa[i,2,k,p]-Pa[i,2,k,p-1])*psi[i,2,k,t,p,l] for p in 2:pow )+Pa[i,2,k,1]*psi[i,2,k,t,1,l]  )*dist*860)/(Masse[k,i]*cpo*r+S[k,i]*dist*K*r)) 
       @constraint(mod,[ i in 1:n, k in 1:Nbp,l in 1:u; JobA[i,2,k,l]>0], y[i,k,tt,2,l]==ti[k,i,l])
      @constraint(mod,[ i in 1:n,k in 1:Nbp,l in 1:u ; JobA[i,2,k,l]>0], temp_min<=y[i,k,ec,2,l]<= temp_max)
      @constraint(mod,[ i in 1:n,k in 1:Nbp,l in 1:u, t in td[k,i,l]:tf[k,i,l]; JobA[i,1,k,l]>0], tmin[k,i,l]<=y[i,k,t,1,l]<= tmax[k,i,l])
      


      @constraint(mod,[i in 1:n,b in 1:B,t in 1:T,l in 1:u],dist*d*Pd[i,b]*(z[b,i,t,l]-1)<=q[b,i,t,l])
      @constraint(mod,[i in 1:n,b in 1:B,t in 1:T,l in 1:u],q[b,i,t,l]<=dist*c*Pc[i,b]*z[b,i,t,l])
      @constraint(mod,[i in 1:n,b in 1:B,t in 1:T-1,l in 1:u],z[b,i,t+1,l]-z[b,i,t,l]<=w[b,i,t,l])
      @constraint(mod,[i in 1:n,b in 1:B,l in 1:u],sum(w[b,i,t,l] for t in 1:T)<=fi)
      
      @constraint(mod,[i in 1:n,b in 1:B,t in 1:T,l in 1:u; gamma[i,b]>0],E[b,i,t,l]==eta*E[b,i,t-1,l]+q[b,i,t,l])
      @constraint(mod,[i in 1:n,b in 1:B],E[b,i,0,1]==xi[i,b])
      @constraint(mod,[i in 1:n,b in 1:B,l in 2:u; u>1 && gamma[i,b]>0],E[b,i,0,l]==E[b,i,T,l-1])
      @constraint(mod,[i in 1:n,t in 1:T,b in 1:B,l in 1:u],E[b,i,t,l]<=gamma[i,b])
      
      @constraint(mod,[i in 1:n,t in 1:T,l in 1:u],sum(f[i,j,t,l] for j in 1:n if i!=j )+Pg[i,t,l]<=P_sous[i])
      @constraint(mod,[i in 1:n,t in 1:T,l in 1:u],sum(f[j,i,t,l] for j in 1:n if i!=j )+Cg[i,t,l]<=P_sous[i])
    
      @constraint(mod,[i in 1:n,l in 1:u],Gain[i,l]==(Pvc*Fsend[i,l]-Pac*Freceive[i,l]+oui*V[i,l]-Pedf*C[i,l]))                       # Gain après
      @constraint(mod,[i in 1:n,l in 1:u],Gain[i,l]>=GainA[i,l]-beta*abs(GainA[i,l]))                                                                # gain après sup a 15% gain avant
    

      @constraint(mod,[i in 1:n,l in 1:u],Fsend[i,l]==dist*sum(f[i,j,t,l] for j in 1:n,t in 1:T if i!=j))
      @constraint(mod,[i in 1:n,l in 1:u],Freceive[i,l]==dist*sum(f[j,i,t,l] for j in 1:n,t in 1:T if i!=j))
      @constraint(mod,[i in 1:n,l in 1:u],V[i,l]==dist*sum(Pg[i,t,l] for t in 1:T ))
      @constraint(mod,[i in 1:n,l in 1:u],C[i,l]==dist*sum(Cg[i,t,l] for t in 1:T ))
      
    optimize!(mod)    
    return round(JuMP.objective_value(mod),digits=2),round(MOI.get(mod, MOI.SolveTimeSec()), digits=2),round(JuMP.objective_bound(mod),digits=2),round(MOI.get(mod, MOI.RelativeGap()),digits=4),1        
end
 
 
 
function appel1(jour::Int,n::Int,Tlimit::Number)
          include("Data/n"*string(n)*"_t48.jl")
          GainA,xa=pasEchange(200.0)
          obj,cpu,Bb,gap,fin=AvecEchange(GainA,jour,200.0)  
          gap*=100

          open("Sol_MILPFortz.txt", "w") do io
              if fin==1
                     if cpu>Tli
                            println(io, "$n & $obj & $gap & $Bb & \$t_l\$")
                     else
                            println(io, "$n & $obj & $gap& $Bb & $cpu")
                     end
              else
                            println(io, "$n & *** & *** & $Bb & \$t_l\$")
             end
          end
end

function appel()
          open("Sol_MILPFortz.txt", "w") do io
          	println(io,"MILP with SOS variables with the day index")
          end
          for n in [7]
                 appel1(1,n,200.0)
         end
         
end
appel()
