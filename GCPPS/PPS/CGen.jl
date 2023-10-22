function CGen(maxIter::Int,GainA::Vector{Float64},xa::Array{Float64},R::Array{Float64},TliS::Number,Tli::Number)                                                                    # Fonction de Génération de colonne 
    P= Array{Float64}(undef, (n,JA,Nbp,T,maxIter))
    a= Array{Float64}(undef, (n,T))
    b= Array{Float64}(undef, (n,JA,Nbp))
    P[:,:,:,:,:].=100
    iter=1
    for i in 1:n
        for j in 1:JA
            for k in 1:Nbp   
                    for t in 1:T
                        P[i,j,k,t,1]=sum(Pa[i,j,k,p]*xa[j,i,k,t,p,1] for p in 1:pow)
                    end
            end 
        end
    end
    C=0
    open("VariationObj_"*string(n)*"_CG_PPs.txt","w") do io
        println(io,"Solution variation file\n")
    end
    
    (mod,alpha,bet,f,cp2)=MasterProb(P,iter,GainA,0,R,Tli)            # Calculate the initial dual values

    C+=cp2
    E=-1000
    
    cout= Array{Float64}(undef, (n,JA,Nbp))                               # reduced costs

    while iter<=maxIter-1 && C<2000
        open("VariationObj_"*string(n)*"_CG_PPs.txt","a") do io
            println(io,"*************************** L'objectif est ",f," à l'itération ",iter,"*************************")
            println(io,"E est ", E)
        end

        if abs(E)<=0.01 && iter>2               # Condition d'arrêt Obj==0
            println("Optimal solution under condition |E|<= epsilon\n")
            return f,P,iter,C
        end
        
        iter+=1
        for i in 1:n
            for t in 1:T
                a[i,t]=dual(alpha[i,t])
            end
        end


        for i in 1:n
            for j in 1:JA
               for k in 1:Nbp
                    if R[i,j,k]>0
                        b[i,j,k]=dual(bet[i,j,k])
                    end
               end
            end
        end
        
        E=0
        for i in 1:n
            for j in 1:JA
                for k in 1:Nbp
                   if R[i,j,k]>0
                        (sch,cout[i,j,k],obj,cp1)=subProb(a[i,:],b[i,j,k],i,j,k,0,TliS) 
                        C=C+cp1
                        if cout[i,j,k]<=0 
                            P[i,j,k,:,iter]=sch[:] 
                            E=E+cout[i,j,k]
                        end
                    end 
                end
            end
        end
        (mod,alpha,bet,f,cp2)=MasterProb(P,iter,GainA,0,R,Tli)  
       
        C=C+cp2
   
        if C>2000
           @info "**********************CPU max atteint*******************************"
           open("VariationObj_"*string(n)*"_CG_PPs.txt","a") do io
                println(io,"E est ", E)
           end
           return f,P,iter,C
        end
    end #while
    @info "La solution optimale n'a pas été obtenue en ",maxIter,"itérations"
    return f,P,iter,C
end
