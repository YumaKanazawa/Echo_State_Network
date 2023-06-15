include("function.jl")
using Plots

#=============教師データ==============#
# max=1.01
# min=0.99

# u=rand(Uniform(min,max),DataLength,K)#一様乱数を取得する
# n=2
# y=NARMA(n,u)

y::Array{Float64,1}=tryparse.(Float64, read_file("../Basic ESN/data/MackeyGlass_t17.txt"))
# y=u

# u=[sin(0.1*i) for i in 1:DataLength]
# y::Array{Float64,1}=u

#各データのスケーリング
# y=y.-(mean(y)-1.0)
#===================================#


#==================各種パラメータ============#
DataLength::Int64=(10000)#length(y)#データの総数
Tr::Int64=4000#学習用データの長さ
T0::Int64=100#過渡時間

Learn_time::Int64=T0+Tr#学習に利用する時間
test_time::Int64=600#テストに利用するデータの長さ

# 各種定数の定義
K::Int64=1#入力のデータ次元
L::Int64=1#出力のデータ次元
N::Int64=600#内部ユニットの個数 Nが小さいと，疎行列の作成ができなくなる(非ゼロの割合で決めているので)
esn=ESN_param(K,N,L)

W_in::Array{Float64,2}=esn.W_in

W_res::Array{Float64,2}=esn.W_rec

W_back::Array{Float64,1}=esn.W_back_V
#===============================#

# Lolenz=RK_N(Lolenz_,Δt,DataLength,[1,1,1])
# y=[]
# for i in 1:DataLength
#     push!(y,Lolenz[i][1])
# end

rn::Array{Float64,1}=zeros(N)#ある時刻でのリザバーベクトル

# γ::Float16=10.0#大きくすれば学習誤差も下がる
# σ::Float16=1.0#分岐点

# println("γ=",γ,",σ=",σ)

X_ESN::Array{Float64,2}=zeros(DataLength,N)#中間層用の格納行列

#=================この上までがリザバーの基本的な設定========================#
# =リザバーの中間層の計算=#
Δt::Float64=17.0/DataLength

function RK_4()
    for n in 1:DataLength
        X_ESN[n:n,:]=rn#時刻nにおける方程式の解

        yn=y[n]

        K1::Array{Float64,1}=γ*Δt*(-rn+tanh.(W_res*rn+σ*W_back*yn))

        K2::Array{Float64,1}=γ*Δt*((-rn+K1/2)+tanh.((W_res*(rn+K1/2))+σ*(W_back*yn)))

        K3::Array{Float64,1}=γ*Δt*(-(rn+K2/2)+tanh.((W_res*(rn+K2/2))+σ*(W_back*yn)))

        K4::Array{Float64,1}=γ*Δt*(-(rn+K3)+tanh.((W_res*(rn+K3))+σ*(W_back*yn)))

        global rn+=(K1+2.0*K2+2.0*K3+K4)/6.0
    end
end

# RK_4()
# Learn_M::Array{Float64,2}=hcat(ones(DataLength),X_ESN)#バイアス項の追加

# λ::Float64=10^-6
# Wout=Learning(Learn_M,y,λ,T0,Learn_time)

# err=NMSE(y[T0:Learn_time],Learn_M[T0:Learn_time,:]*Wout)
# println("err=",err)
# plot(Learn_M[T0:Learn_time,:]*Wout,xlim=(T0+100,T0+200))
# plot!(y[T0:Learn_time,:])

#============================予測値の計算==================================#
r_new=X_ESN[Learn_time:Learn_time,:]'#学習の最終時刻のデータをベクトル化
function RK_4_test()
    for n in Learn_time+1:DataLength
        yn=(hcat(1.0,r_new')*Wout)[1]

        K1=γ*Δt*(-r_new+tanh.(W_res*r_new+W_back*yn))

        K2=γ*Δt*(-(r_new+K1/2)+tanh.(W_res*(r_new+K1/2)+σ*W_back*yn))

        K3=γ*Δt*(-(r_new+K2/2)+tanh.(W_res*(r_new+K2/2)+σ*W_back*yn))

        K4=γ*Δt*(-(r_new+K3)+tanh.(W_res*(r_new+K3)+σ*W_back*yn))

        global r_new+=(K1+2.0*K2+2.0*K3+K4)/6.0

        X_ESN[n:n,:]=r_new#時刻nにおける方程式の解
    end
end

# RK_4_test()

# Pred=hcat(ones(DataLength),X_ESN)
# Pred=X_ESN

# err_NMSE=NMSE(y[Learn_time+1:Learn_time+1+test_time],Pred[Learn_time+1:Learn_time+1+test_time,:]*Wout)
# println("Err=",err_NMSE)


# T=Learn_time
# plot(Pred*Wout,label=("ESN"),xlims=(T-100,T+100),ylims=(minimum(y),maximum(y)))
# plot!(y)
# savefig("macky_glass.png")


#===================分岐解析=====================#
function δ(i::Int64,j::Int64)
    ret=0.0;
    if i==j
        ret=1.0
    else
        ret=0.0
    end
    return ret
end

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

    J::Array{Float64,2}=zeros(N,N)
    s=W_res*rn+σ*W_back*y[DataLength]
    for i in 1:N
        for j in 1:N
            J[i,j]=γ*(δ(i,j)+W_res[i,j]/cosh(s[i]))
        end
    end

    # #全ての固有値が正であるか，負であるか(安定性解析)
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
    how_many_positive::Int64=0
    for i in 1:length(eigJ)
        if eigJ[i]>0.0
            how_many_positive+=1
        end
    end
    println("σ=",Σ,",",how_many_positive)
    global Σ+=0.001
end




