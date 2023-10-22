
function CGen(maxIter::Int,GainA::Vector{Float64},xa::Array{Float64},R::Array{Float64},P::Array{Float64},Tli::Float64,TliS::Float64)
    open("VariationObj_"*string(n)*"_CG_PPoInit.txt","w") do io
    end
    P[:,:,:,:,:].=Inf
    iter=1	
    C=0
    maxTimeCG=2000
    cp2=0
    f=zeros(n,n,T)
    
    for i in 1:n
        for j in 1:JA
            for k in 1:Nbp   
                  for t in 1:T
                               P[i,j,k,t,iter]=sum(Pa[i,j,k,p]*xa[j,i,k,t,p,1] for p in 1:pow)
                end
             end 
        end
    end
    
    (mod,alpha,bet,f,cp2)=MasterProb(P,iter,GainA,0,R,Tli) 

    C=C+cp2
    obj=-1000
    gap=0
   
    while iter<=maxIter-1 && C<maxTimeCG
    
        open("VariationObj_"*string(n)*"_CG_PPoInit.txt","a") do io
            println(io,"*************************** The objective equals ",f," at iteration ",iter,"*************************")
            println(io,"E equals ", obj ," with gap from the best bound ",gap)
        end
     
        if abs(obj)<=0.01               # Condition d'arrêt Obj==0
            println("Optimal solution founded under condition |E|<= epsilon\n")
                 return f,P,iter,C
        end
        
        iter=iter+1

       (sch,E,obj,cp1,gap)=subProb(alpha,bet,0,R,TliS) 

        for i in 1:n
            for j in 1:JA
                for k in Nbp
                    if E[i,j,k]<=0 && R[i,j,k]>0
                                   P[i,j,k,:,iter]=view(sch,i,j,k,:)
                    end
                end
            end
        end
        (mod,alpha,bet,f,cp2)=MasterProb(P,iter,GainA,0,R,Tli)
        C=C+cp1+cp2
    end     #end while
    
    println("The optimal solution was not found after ",maxIter," iterations")
    iter-=1
    return f,P,iter,C

end  #end function

