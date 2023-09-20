### By Mariam SANGARE ###
### November 2021     ###
using JuMP
using CPLEX

Tli=5600
TliS=200
include("pasEchange.jl")                 # Determines the initial solution
include("Maitre.jl")                           # Master problem
include("AfficheGc.jl")                     # print the output and the graphs in latex format
include("subProblemPPo.jl")          # subproblem
include("CGenPPO.jl")                     # column generation based heuristic

function MainPPo() 
   open("Sol_PPsFortz.txt", "w") do io
   		println(io,"Solutions in latex format")
    end
   #for n in [7 28 56 112 224]
  for n in [7]
           include("n"*string(n)*"_t48.txt")
           P= Array{Float64}(undef, (n,JA,Nbp,T,maxIterGC))
           R= Array{Float64}(undef, (n,JA,Nbp))             # R[i,j,k] is greater than one if task j is required in room k at least once over the planning horizon
           R[:,:,:].=sum(JobA[:,:,:,l] for l in 1:u)

           GainA,xa=pasEchange()                       
           f_e,P,iter,cpu1=CGen(maxIterGC,GainA,xa,R,P)
        
           f,Ce,Ve,Gain,tote,inje,exte,Besse,xe,Puise,cpu2,q=MasterProb(P,maxIterGC,GainA,1,R)
           cpuf=cpu1+cpu2
           
          open("VariationObj_"*string(n)*"_CG_PPoInit.txt","a") do io
          		println(io,"\n solution is $f the cpu time : $cpuf\n")
           end
           #aa=affiche(f,Ce,Ve,Gain,tote,inje,exte,Besse,xe,Puise,cpuf,q)         # print the graph in latex format
            
           f=round(f,digits=2)
           open("Sol_PPoInit.txt", "a") do io
                         println(io, "$f & $cpuf & gap")
          end
    end
end
MainPPo()
