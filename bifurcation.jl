include("ESN_setting.jl")

#===================分岐解析=====================#
function Jacobi_real_eigvals(σ::Float64,γ::Float64)
    rn::Array{Float64,1}=zeros(N)
    for n in 1:10^4
        yn=y[n]

        K1::Array{Float64,1}=γ*Δt*(-rn+tanh.(W_res*rn+σ*W_back*yn))

        K2::Array{Float64,1}=γ*Δt*((-rn+K1/2)+tanh.((W_res*(rn+K1/2))+σ*(W_back*yn)))

        K3::Array{Float64,1}=γ*Δt*(-(rn+K2/2)+tanh.((W_res*(rn+K2/2))+σ*(W_back*yn)))

        K4::Array{Float64,1}=γ*Δt*(-(rn+K3)+tanh.((W_res*(rn+K3))+σ*(W_back*yn)))
        
        rn+=(K1+2.0*K2+2.0*K3+K4)/6.0
    end

    J=Jacobi_matrix(γ,σ)
    #全ての固有値が正であるか，負であるか(安定性解析)
    # eigJ::Array{Float64,1}=real.(eigvals(J))
    # how_many_positive::Int64=0
    # for i in 1:length(eigJ)
    #     if eigJ[i]>0.0
    #         global how_many_positive+=1
    #     end
    # end

    return real.(eigvals(J))
end

Σ::Float64=0.0
while Σ<10.0
    #全ての固有値が正であるか，負であるか(安定性解析)
    eigJ::Array{Float64,1}=Jacobi_real_eigvals(Σ,10.0)
    local how_many_positive::Int64=0
    for i in 1:length(eigJ)
        if eigJ[i]>0.0
            how_many_positive+=1
        end
    end
    println("σ=",Σ,",",how_many_positive)
    global Σ+=0.001
end