function affiche(f,Ce,Ve,Gain,tote,inje,exte,Besse,xe,Puise,cpuf,q)
 


     open("VariationObj_"*string(n)*"_CG_plusieurs jours_pricing une fois.txt","a") do io
          println(io,"\n L'objectif après possiblité d'échanges ", f,"kWh \n")
          println(io," Les gains après possiblité d'échanges ", Gain,"kWh\n")
          println(io," Les soutirages du réseau après possiblité d'échanges ", Ce,"kWh\n")
          #println(io," Les injections dans le réseau après possiblité d'échanges ", Ve ,"kWh\n")
          println(io," La demande totale après possiblité d'échanges ", tote,"kWh\n")
          println(io," Les injections périodiques après possiblité d'échanges ", inje,"kWh\n")
          println(io," Les extractions périodiques après possiblité d'échanges ", exte,"kWh\n")
          println(io," Les niveaux périodiques des batteries après possiblité d'échanges ", Besse,"kWh\n")
          println(io," Le temps CPU est $cpuf")
    end


   bat=0

    for i in 1:n
                 bat=bat+xi[i]
    end
    production= Array{Float64}(undef, (T))
    for t in 1:T
          production[t]=sum(prod1[i,t] for i in 1:n)
    end

    s=""
    for t in 1:T
                      s=s*"("*string(t)*","*string(production[t])*")"
    end

    open("VariationObj_"*string(n)*"_CG_plusieurs jours_pricing une fois.txt","a") do io
          println(io," \n******************APRÈS************************\n")
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
                legend entries={Production(kW),Demande totale(kWh), Soutirée(kWh),Injectée(kWh),BESS(kWh)},
                legend style={at={(0,1)},anchor=north west}]")
          println(io,"\\addplot[green, very thick,mark=*] coordinates{",s,"};")
    end

    s=" "
for t in 1:T
    s=s*"("*string(t)*","*string(tote[t])*")"
end
          open("VariationObj_"*string(n)*"_CG_plusieurs jours_pricing une fois.txt","a") do io
        println(io,"\\addplot[orange, very thick,mark=+] coordinates{",s,"};")
end

s=" "
for t in 1:T
    s=s*"("*string(t)*","*string(exte[t])*")"
end
          open("VariationObj_"*string(n)*"_CG_plusieurs jours_pricing une fois.txt","a") do io
          println(io,"\\addplot[blue, very thick,mark=o] coordinates{",s,"};")
end



s=" "
for t in 1:T
    s=s*"("*string(t)*","*string(inje[t])*")"
end
          open("VariationObj_"*string(n)*"_CG_plusieurs jours_pricing une fois.txt","a") do io
          println(io,"\\addplot[red, very thick,mark=x] coordinates{",s,"};")
end

#=
s="(0,"*string(bat)*")"
for t in 1:T
    s=s*"("*string(t)*","*string(Besse[t])*")"
end
=#
BESS_inj= Array{Float64}(undef, (T))
       s=" "
   for t in 1:T
          BESS_inj[t]=sum(q[b,i,t] for i in 1:n,b in 1:B if q[b,i,t]>=0)
          s=s*"("*string(t)*","*string(BESS_inj[t])*")"

   end

          open("VariationObj_"*string(n)*"_CG_plusieurs jours_pricing une fois.txt","a") do io
          println(io,"\\addplot[brown, very thick,mark=|] coordinates{",s,"};")
          println(io,"  \\end{axis}
          \\end{tikzpicture}    
          \\caption{ Avec échanges entre les membres.}
            \\end{figure}")
    end

end



