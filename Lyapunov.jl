include("ESN_setting.jl")

#=リアプノフ指数の計算=#
function Lyapnov(γ,σ)
    Linear_of_govern_eq(t,x)=Jacobi_matrix(γ,σ)*x#線形化方程式
    #微小な摂動は線形化方程式を満たす.
    δ0=ones(N)
    iterate::Int16=1000#反復回数
    solution=RK_N(Linear_of_govern_eq,Δt,iterate,δ0)

    ryapnov::Array{Float32,1}=zeros(N)
    
    for k in 1:iterate-1
        ratio::Array{Float32,1}=abs.(solution[k+1])./abs.(solution[k])
        ryapnov+=(1/k).*log.(ratio)
    end
    return ryapnov
end

Lyapnov(10.0,0.01)