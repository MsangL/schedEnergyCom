### By Mariam SANGARE ###
### November 2021     ###

using JuMP, CPLEX

function main()
  include("pasEchange.jl")
  include("Maitre.jl") 
  include("CGen.jl")
  include("subProblem.jl")
  #include("AfficheGc.jl")
  Tli=3600.0
  TliS=50.0
  open("Sol_PPs.jl", "w") do io
      println(io,"This file contains the solutions\n")
  end
  # for n in [112 224]
   for n in [7]
           include("Data/n"*string(n)*"_t48.jl")
           R= Array{Float64}(undef, (n,JA,Nbp))             # sera supérieure à 0 si la tâche est demandée au moins 1 fois
           for i in 1:n
                  for j in 1:JA
                      for k in 1:Nbp
                           R[i,j,k]=sum(JobA[i,j,k,l] for l in 1:u)
                     end
                 end
            end
           GainA,xa=pasEchange(TliS)                       
           f_e,P,iter,cpu1=CGen(maxIterGC,GainA,xa,R,TliS,Tli)
           f,Ce,Ve,Gain,tote,inje,exte,Besse,xe,Puise,cpu2,q=MasterProb(P,maxIterGC,GainA,1,R,Tli)
           cpuf=cpu1+cpu2
           
          open("VariationObj_"*string(n)*"_CG_PPs.jl","a") do io
            println(io,"\n La solution obtenue est $f et le temps est : $cpuf et iter finale est $iter\n")
          end
         #affiche(f,Ce,Ve,Gain,tote,inje,exte,Besse,xe,Puise,cpuf,q)
           f=round(f,digits=2)
          open("Sol_PPs.jl", "a") do io
            println(io, "$f & $cpuf & gap")
          end
    end
end
main()


