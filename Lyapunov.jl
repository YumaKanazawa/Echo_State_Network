include("ESN_setting.jl")


#=リアプノフ指数の計算=#
function Lyapnov(γ,σ)
    println("γ=",γ,",σ=",σ)
    Δt=0.01#時間ステップ
    T=1000#反復回数
    #時刻nにおけるJacobi行列
    xn::Array{Float64,1}=zeros(N)
    for n in 1:T

        if n==1
            yn=0.0
        else
            yn=y[n-1]
        end

        tn::Float64=n*Δt

        K1::Array{Float64,1}=γ*Δt*(-rn+tanh.(W_res*rn+σ*W_back*yn))

        K2::Array{Float64,1}=γ*Δt*((-rn+K1/2)+tanh.((W_res*(rn+K1/2))+σ*(W_back*yn)))

        K3::Array{Float64,1}=γ*Δt*(-(rn+K2/2)+tanh.((W_res*(rn+K2/2))+σ*(W_back*yn)))

        K4::Array{Float64,1}=γ*Δt*(-(rn+K3)+tanh.((W_res*(rn+K3))+σ*(W_back*yn)))

        xn+=(K1+2.0*K2+2.0*K3+K4)/6.0
    end

    Jn::Array{Float64,2}=Jacobi_matrix(γ,σ,T,xn)#時刻nにおけるJacobi行列
    eigJ=eigvals(Jn'*Jn).^(1/(2*T))#Jnの固有値

    Lyapnov_diameter=(1/T).*log.(eigJ)
    return maximum(Lyapnov_diameter)
end


#=リアプノフ指数の計算=#
function Lyapnov_lorentz(r)
    println("r=",r)
    Δt=0.01#時間ステップ
    T=10000#反復回数
    #時刻nにおけるJacobi行列
    xn::Array{Float64,1}=[10.0,10.0,10.0]
    # δn::Array{Float64,1}=zeros(3)#微小な摂動

    ret::Array{Float64,1}=zeros(3)

    for n in 1:T

        tn::Float64=n*Δt

        K1::Array{Float64,1}=Δt*(Lolenz_(tn,xn,r))

        K2::Array{Float64,1}=Δt*(Lolenz_(tn+Δt/2,xn+K1/2,r))

        K3::Array{Float64,1}=Δt*(Lolenz_(tn+Δt/2,xn+K2/2,r))

        K4::Array{Float64,1}=Δt*(Lolenz_(tn+Δt,xn+K3,r))


        #=リアプノフ指数の計算=#
        #=微小な摂動δ(t)はJacobianを用いてδ'(t)=J(t)δ(t)に従う=#
        Jn::Array{Float64,2}=zeros(3,3)#ローレンツ方程式のやこび行列

        Jn[1,1]=-10.0
        Jn[1,2]=10.0
        Jn[1,3]=0.0
        Jn[2,1]=r-xn[3]
        Jn[2,2]=-1.0
        Jn[2,3]=-xn[1]
        Jn[3,1]=xn[2]
        Jn[3,2]=xn[1]
        Jn[3,3]=-(8/3)

        # K1_l=Δt*Jn*δn
        # K2_l=Δt*Jn*(δn+K1_l/2)
        # K3_l=Δt*Jn*(δn+K2_l/2)
        # K4_l=Δt*Jn*(δn+K3_l)

        # δnew=δn+(K1_l+2.0*K2_l+2.0*K3_l+K4_l)/6.0
        # logabs.(δnew)./abs.(δ)
        
        # δn=δnew

        eigJ=(eigvals(Jn'*Jn)).^(1/(2*tn))
        ret=log.(eigJ)/tn
        #==================#

        xn+=(K1+2.0*K2+2.0*K3+K4)/6.0
    end
    return maximum(ret)
end

Λ=[r for r in 18.0:0.01:26.0]
plotting=[Lyapnov_lorentz(r) for r in 18.0:0.01:26.0]

plot(Λ,plotting)