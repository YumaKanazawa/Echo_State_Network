include("ESN_setting.jl")


#=固定点の探索=#
function Newton_N(F,J)
    d=N
    x::Array{Float64,1}=[-0.001 for i in 1:d]#初期値の設定
    i=0
    while true
        x_new=x-inv(J(x))*F(x)
        if norm(x-x_new,2.0)<10^-10 || i>=10^10
            break
        end
        x=x_new
        i+=1
    end
    return x
end


function fix_point(γ::Float64,σ::Float64,n::Int64)
    F(x)=γ*(-x+tanh.(W_res*x+σ*W_back*y[n]))
    J(x)=Jacobi_matrix(γ,σ,n,x)

    # F(x)=[σ*x[1]-γ*x[2],σ*x[2]+γ*x[1]]
    # J(x)=[σ -γ;γ σ]

    fix=Newton_N(F,J)#固定点

    J_fix=J(fix)
    real_eigval=real.(eigvals(J_fix))

    return real_eigval
end



