using JuMP, CPLEX
u=3
include("pasEchange.jl")
include("MB_js.jl")
function solveOutput(jour::Int,Tli::Number,Tlimit)

    include("Data/jour"*string(jour)*".jl")
    Tli=200.0
     C0,GainA,E0,y0,gap1,inj,ext,Bess,tot,q0=PasEchange(jour,Tli)
     somme=0
     for i in 1:n
         somme=somme+C0[i]
     end

     s=""
     for b in 1:B
          for i in 1:n
              s=s*" "*string(E0[b,i,T])*";"
          end
     end
     
     open("Data/jour"*string(jour+1)*".jl","a") do io
            println(io,"xi0=["*chop(s)*"]")
     end
     
     s=""
     for k in 1:Nbp
          for i in 1:n
              s=s*" "*string(y0[i,k,T,1])
          end
          s=s*";"
     end
          
     
     open("Data/jour"*string(jour+1)*".jl","a") do io
            println(io,"y0=["*chop(s)*"]")
     end
     
     s=""
     for k in 1:Nbp
          for i in 1:n
              s=s*" "*string(y0[i,k,T,2])
          end
          s=s*";"
     end
          
     
     open("Data/jour"*string(jour+1)*".jl","a") do io
            println(io,"te0=["*chop(s)*"]")
     end
    #*********************************************************************
     cpu,C,Gain,E,y,gap,inje,exte,Besse,tote,q,xa=AvecEchanges(GainA,jour,Tlimit)


     somme=0
     for i in 1:n
         somme=somme+C[i]
     end

     production=zeros(T)
     for t in 1:T
          production[t]=sum(prod1[i,t] for i in 1:n)
     end

     s=""
     for t in 1:T
           s=s*"("*string((jour-1)*48+t)*","*string(production[t])*")"
     end

     open("Solutions_n"*string(n)*"_jour"*string(jour)*"SansInd.jl","a") do io
               println(io," ******************Figures before members exchanges************************\n The objective is $C0 kWh at day $jour \n")
               println(io,"
               \\begin{figure}[!h]
               \\centering
               \\begin{tikzpicture}
               \\begin{axis}[grid= major ,
               width=0.6\\textwidth ,
               xlabel = {\$Temps\$ (\$\\frac{1}{2}h\$)} ,
               ylabel = {} ,
               xmin = 0, xmax = 48,
               ymin = 0, ymax = 30,
               legend entries={Production(kW),Total demand (kWh),Injected(kWh),BESS(kWh)},
               legend style={at={(0,1)},anchor=north west}]")
               println(io,"\\addplot coordinates{",s,"};")
    end

    s=" "
    for t in 1:T
          s=s*"("*string((jour-1)*48+t)*","*string(tot[t])*")"
    end
     open("Solutions_n"*string(n)*"_jour"*string(jour)*"SansInd.jl","a") do io
          println(io,"\\addplot coordinates{",s,"};")
    end

    s=" "
    for t in 1:T
          s=s*"("*string((jour-1)*48+t)*","*string(ext[t])*")"
    end
     open("Solutions_n"*string(n)*"_jour"*string(jour)*"SansInd.jl","a") do io
          println(io,"\\addplot coordinates{",s,"};")
    end

    s=" "
    for t in 1:T
          s=s*"("*string((jour-1)*48+t)*","*string(inj[t])*")"
    end

     open("Solutions_n"*string(n)*"_jour"*string(jour)*"SansInd.jl","a") do io
             println(io,"\\addplot coordinates{",s,"};")
    end

    s=""
    for t in 1:T
    somme=0
        for i in 1:n
              for b in 1:B
                  if (q0[b,i,t]>0)
                       somme=somme+q0[b,i,t]
                  end
              end
        end  
    s=s*"("*string((jour-1)*48+t)*","*string(somme)*")"
    end
     open("Solutions_n"*string(n)*"_jour"*string(jour)*"SansInd.jl","a") do io
             println(io," \\addplot+[color=green] coordinates{",s,"};")
             println(io,"  \\end{axis}
             \\end{tikzpicture}    
             \\caption{Sans échanges entre les membres.}
             \\end{figure}")
    end



     
     s=""
     for b in 1:B
          for i in 1:n
              s=s*" "*string(E[b,i,T])*";"
          end
     end
     
     open("Data/jour"*string(jour+1)*".jl","a") do io
            println("##Batteries et températures avec échanges")
            println(io,"xi=["*chop(s)*"]")
     end
     
     s=""
     for k in 1:Nbp
          for i in 1:n
              s=s*" "*string(y[i,k,T,1])
          end
          s=s*";"
     end
          
     
     open("Data/jour"*string(jour+1)*".jl","a") do io
            println(io,"y00=["*chop(s)*"]")
     end
     
     s=""
     for k in 1:Nbp
          for i in 1:n
              s=s*" "*string(y[i,k,T,2])
          end
          s=s*";"
     end
          
     
     open("Data/jour"*string(jour+1)*".jl","a") do io
            println(io,"te=["*chop(s)*"]")
     end
     
     
     production= Array{Float64}(undef, (T))
     for t in 1:T
          production[t]=sum(prod1[i,t] for i in 1:n)
     end

     s=""
     for t in 1:T
           s=s*"("*string((jour-1)*48+t)*","*string(production[t])*")"
     end

     open("Solutions_n"*string(n)*"_jour"*string(jour)*"SansInd.jl","a") do io
               println(io," ****************** When members exchange their surplus the objective is $C kWh at day $jour ************************\n")
               println(io,"
               \\begin{figure}[!h]
               \\centering
               \\begin{tikzpicture}
               \\begin{axis}[grid= major ,
               width=0.6\\textwidth ,
               xlabel = {\$Temps\$ (\$\\frac{1}{2}h\$)} ,
               ylabel = {} ,
               xmin = 0, xmax = 48,
               ymin = 0, ymax = 30,
               legend entries={Production(kW),Total demand(kWh), Withdrawn (kWh),Injected(kWh),BESS(kWh)},
               legend style={at={(0,1)},anchor=north west}]")
               println(io,"\\addplot coordinates{",s,"};")
    end

    s=" "
    for t in 1:T
          s=s*"("*string((jour-1)*48+t)*","*string(tote[t])*")"
    end
     open("Solutions_n"*string(n)*"_jour"*string(jour)*"SansInd.jl","a") do io
          println(io,"\\addplot coordinates{",s,"};")
    end

    s=" "
    for t in 1:T
          s=s*"("*string((jour-1)*48+t)*","*string(exte[t])*")"
    end
     open("Solutions_n"*string(n)*"_jour"*string(jour)*"SansInd.jl","a") do io
          println(io,"\\addplot coordinates{",s,"};")
    end

    s=" "
    for t in 1:T
          s=s*"("*string((jour-1)*48+t)*","*string(inje[t])*")"
    end

     open("Solutions_n"*string(n)*"_jour"*string(jour)*"SansInd.jl","a") do io
             println(io,"\\addplot coordinates{",s,"};")
    end

    s=""
    for t in 1:T
    somme=0
        for i in 1:n
              for b in 1:B
                  if (q[b,i,t]>0)
                       somme=somme+q[b,i,t]
                  end
              end
        end  
    s=s*"("*string((jour-1)*48+t)*","*string(somme)*")"
    end
     open("Solutions_n"*string(n)*"_jour"*string(jour)*"SansInd.jl","a") do io
             println(io," \\addplot+[color=green] coordinates{",s,"};")
             println(io,"  \\end{axis}
             \\end{tikzpicture}    
             \\caption{Avec échanges entre les membres.}
             \\end{figure}")
    end
    obj=sum(C[i] for i in 1:n)
    return obj,cpu

end

function main()
    C=0
    cpu=0

    for j in 1:u
        c1,c2=solveOutput(j,50,50)
        cpu=cpu+c2
        C=C+c1
    end
    open("Solution_"*string(u)*"days.jl") do io
        println(io,"obj $C and cpu : $cpu")
    end

end

main()