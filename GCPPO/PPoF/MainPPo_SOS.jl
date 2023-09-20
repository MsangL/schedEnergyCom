### By Mariam SANGARE ###
### November 2021     ###
using JuMP
using CPLEX

Tli=5600
TliS=200
include("pasEchange.jl")
include("Maitre.jl")
include("AfficheGc.jl")
include("subProblemPPO_SOS.jl")
include("CGPPO_SOS.jl")

lim=0.6


function appel() 
   open("Sol_PPoFortz.txt", "w") do io
   end
# for n in [224]
   for n in [7]
           include("n"*string(n)*"_t48.txt")
           P= Array{Float64}(undef, (n,JA,Nbp,T,maxIterGC))
           R= Array{Float64}(undef, (n,JA,Nbp))            
           R[:,:,:].=sum(JobA[:,:,:,l] for l in 1:u)
           
           GainA,xa=pasEchange()                       # Determine the initial solution
           f_e,P,iter,cpu1=CGen(maxIterGC,GainA,xa,R,P)
           f,Ce,Ve,Gain,tote,inje,exte,Besse,xe,Puise,cpu2,q=MasterProb(P,maxIterGC,GainA,1,R)
           cpuf=cpu1+cpu2
           
          open("VariationObj_"*string(n)*"_CG_PPoFortz.txt","a") do io
                   println(io,"\n La solution obtenue est $f et le temps est : $cpuf et la dernière itération était à : $iter\n")
           end
           #aa=affiche(f,Ce,Ve,Gain,tote,inje,exte,Besse,xe,Puise,cpuf,q)
            
           f=round(f,digits=2)
           open("Sol_PPoFortz.txt", "a") do io
                         println(io, "$f & $cpuf & gap")
          end
    end
end
appel()
