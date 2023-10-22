### By Mariam SANGARE ###
### November 2021     ###
using JuMP,CPLEX

#Tli : is the master problem's solving time limit
#TliS : is the subproblem's solving time limit


function MainPPo(Tli::Number,TliS::Number) 
    include("pasEchange.jl")                 # Determines the initial solution
    include("Maitre.jl")                     # Master problem
    include("AfficheGc.jl")                  # print the output and the graphs in latex format
    include("subProblemPPo.jl")              # subproblem
    include("CGenPPO.jl")                    # column generation based heuristic

    TimeLimitInitialSolution=TliS
    open("Sol_PPoInit.txt", "w") do io
   		println(io,"This file contains the outputs in latex format\n")
    end
   #for n in [7 28 56 112 224]
  for n in [7]
        include("Data/n"*string(n)*"_t48.jl")
        P= Array{Float64}(undef, (n,JA,Nbp,T,maxIterGC))
        R= Array{Float64}(undef, (n,JA,Nbp))             # R[i,j,k] is greater than one if task j is required in room k at least once over the planning horizon
        xa=zeros(JA,n,Nbp,T,pow)
        GainA=zeros(n)
        view(R,:,:,:).=sum(view(JobA,:,:,:,l) for l in 1:u)
        cpu1=0
        cpuf=0
        GainA,xa=pasEchange(TliS)                          # Determines the initial solution    
        f,P,iter,cpu1=CGen(maxIterGC,GainA,xa,R,P,Tli,TliS)
        cpuf=cpuf+cpu1
        f,Ce,Ve,Gain,tote,inje,exte,Besse,xe,Puise,cpu1,q=MasterProb(P,maxIterGC,GainA,1,R,Tli) # solve master problem with integrality constrains
        cpuf=cpuf+cpu1
        open("VariationObj_"*string(n)*"_CG_PPoInit.txt","a") do io
          	println(io,"\n solution is $f the cpu time : $cpuf\n")
        end
        #affiche(f,Ce,Ve,Gain,tote,inje,exte,Besse,xe,Puise,cpuf,q)         # print the graphs in latex format 
            
        open("Sol_PPoInit.txt", "a") do io
            println(io,  round(f,digits=2),"& $cpuf & gap")
        end
    end
end
MainPPo(3600,200)

