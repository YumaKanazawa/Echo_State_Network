include("ESN_setting.jl")

τ=20

λ::Float64=10^-6

#=================この上までがリザバーの基本的な設定========================#
function training(γ,σ)
    local X_ESN::Array{Float64,2}=zeros(DataLength,N)#中間層用の格納行列
    local Ans::Array{Float64,1}=zeros(DataLength)#正解データの格納
    rn::Array{Float64,1}=zeros(N)#リザバーの初期値
    
    k=1#格納する位置の参照
    n=1
    while k<=Learn_time

        if n==k*τ#Δtのτ倍(例:Δt=0.001の時，0.02(=20*Δt)ごとにデータを格納していく)
            k+=1
        end

        X_ESN[n:n,:]=rn#時刻nにおける方程式の解
        Ans[n]=y[n]

        #時間に応じて予測するデータの反映の仕方を決める
        #学習用のデータの反映
        # if n==1
        #     yn=0.0
        # else
        #     yn=y[n-1]
        # end

        yn=y[n]

        println("L:k=",k,",n=",n,",後:",Learn_time*τ-n)

        K1::Array{Float64,1}=γ*Δt*(-rn+tanh.(W_res*rn+σ*W_back*yn))

        K2::Array{Float64,1}=γ*Δt*(-(rn+K1/2)+tanh.((W_res*(rn+K1/2))+σ*(W_back*yn)))

        K3::Array{Float64,1}=γ*Δt*(-(rn+K2/2)+tanh.((W_res*(rn+K2/2))+σ*(W_back*yn)))

        K4::Array{Float64,1}=γ*Δt*(-(rn+K3)+tanh.((W_res*(rn+K3))+σ*(W_back*yn)))

        rn+=(K1+2.0*K2+2.0*K3+K4)/6.0
        n+=1
    end

    Learn_M::Array{Float64,2}=hcat(ones(DataLength),X_ESN,X_ESN.^2)

    Wout::Array{Float64,1}=Learning(Learn_M,Ans,λ,T0,Learn_time)

    # k=1
    while n<=DataLength
        #予測のデータ反映について
        # if n==Learn_time*τ+1
        #     yn=y[n]
        # else
        #     yn::Float64=hcat(1.0,rn',rn'.^2)⋅Wout
        # end

        yn::Float64=hcat(1.0,rn',rn'.^2)⋅Wout

        X_ESN[n:n,:]=rn#時刻nにおける方程式の解

        # if n==k*τ
        #     # Ans[k]=yn
        #     k+=1
        # end

        println("T:n=",n,",Err=",abs(y[n]-yn))

        K1::Array{Float64,1}=γ*Δt*(-rn+tanh.(W_res*rn+σ*W_back*yn))

        K2::Array{Float64,1}=γ*Δt*(-(rn+K1/2)+tanh.((W_res*(rn+K1/2))+σ*(W_back*yn)))

        K3::Array{Float64,1}=γ*Δt*(-(rn+K2/2)+tanh.((W_res*(rn+K2/2))+σ*(W_back*yn)))

        K4::Array{Float64,1}=γ*Δt*(-(rn+K3)+tanh.((W_res*(rn+K3))+σ*(W_back*yn)))

        rn+=(K1+2.0*K2+2.0*K3+K4)/6.0
        n+=1
    end

    T=Learn_time*τ
    Pred=hcat(ones(DataLength),X_ESN,X_ESN.^2)

    # err_p=NMSE(y,pred)
    err=NMSE(y[Learn_time*τ:DataLength],(Pred*Wout)[Learn_time*τ:DataLength])

    println("LearnErr=",NMSE(y[T0*τ:Learn_time*τ],(Pred*Wout)[T0*τ:Learn_time*τ]))
    println("LearnErr=",NMSE(y[Learn_time*τ:DataLength],(Pred*Wout)[Learn_time*τ:DataLength]))
    plot(Pred*Wout,xlim=(T,DataLength),ylim=(minimum(y),maximum(y)),title=("τ="*string(τ)*",σ="*string(σ)*",γ="*string(γ)*",Err="*string(err)))
    plot!(y)
    savefig("image/σ="*string(σ)*".png")

    # return (Pred*Wout)[T:DataLength]
end


function ESP(γ,σ,rn::Array{Float64,1},r1::Array{Float64,1})
    println("γ=",γ,",σ=",σ)
    n=1
    T=5000
    while n<=T
        #時間に応じて予測するデータの反映の仕方を決める
        #学習用のデータの反映
        if n==1
            yn=0.0
        else
            yn=y[n-1]
        end

        # yn=y[n]

        #================初期値rnにおける数値解===================#
        K1::Array{Float64,1}=γ*Δt*(-rn+tanh.(W_res*rn+σ*W_back*yn))

        K2::Array{Float64,1}=γ*Δt*((-rn+K1/2)+tanh.((W_res*(rn+K1/2))+σ*(W_back*yn)))

        K3::Array{Float64,1}=γ*Δt*(-(rn+K2/2)+tanh.((W_res*(rn+K2/2))+σ*(W_back*yn)))

        K4::Array{Float64,1}=γ*Δt*(-(rn+K3)+tanh.((W_res*(rn+K3))+σ*(W_back*yn)))

        rn+=(K1+2.0*K2+2.0*K3+K4)/6.0
        #======================================================#

        #================初期値r1における数値解===================#
        K1_1::Array{Float64,1}=γ*Δt*(-r1+tanh.(W_res*r1+σ*W_back*yn))

        K2_1::Array{Float64,1}=γ*Δt*((-r1+K1_1/2)+tanh.((W_res*(r1+K1_1/2))+σ*(W_back*yn)))

        K3_1::Array{Float64,1}=γ*Δt*(-(r1+K2_1/2)+tanh.((W_res*(r1+K2_1/2))+σ*(W_back*yn)))

        K4_1::Array{Float64,1}=γ*Δt*(-(r1+K3_1)+tanh.((W_res*(r1+K3_1))+σ*(W_back*yn)))

        r1+=(K1_1+2.0*K2_1+2.0*K3_1+K4_1)/6.0
        #======================================================#

        println("n=",n,",Err=",NMSE(rn,r1))#時刻nにおける異なる初期値の誤差

        n+=1
    end 
end


#=周期解の確認=#
training(10.0,0.012)


# init::Array{Float64,1}=ones(N)#初期値
# init1::Array{Float64,1}=[10.0 for i in 1:N]#初期値からの微妙なずれ

# ESP(10.0,0.012,init,init1)
