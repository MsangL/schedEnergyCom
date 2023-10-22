using JuMP, CPLEX


function appel()
	include("indice.jl")
	include("pasEchange.jl")
	Tli=200.0
	Tlimit=5600.0
          open("Sol_MILPInd.jl", "w") do io
          		println(io,"Solutions\n")
          end
          
          for n in [232425]
              #for n in [23 2324  232425]
                     include("Data/Ind"*string(n)*".jl")      
                GainA,xa=pasEchange(Tli)
                obj,cpu,Bb,gap,DAYC=AvecEchange(GainA,Tlimit)  
                gap*=100
                
                open("Sol_MILPInd.jl", "a") do io
                         if cpu>Tli       
                         println(io, "$n & $obj & $gap & $Bb & \$t_l\$ the consumtion per day is $DAYC")
                         else
                                println(io, "$n & $obj & $gap& $Bb & $cpu the consumption per day is $DAYC")
                         end
                 end
          end
end

appel()


