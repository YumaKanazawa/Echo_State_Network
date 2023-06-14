include("function.jl")
using Plots

#=============教師データ==============#
# max=1.01
# min=0.99

# u=rand(Uniform(min,max),Length,K)#一様乱数を取得する
# n=2
# y=NARMA(n,u)

# y=tryparse.(Float64, read_file("Basic ESN/data/MackeyGlass_t17.txt"))
# y=u

# u=[sin(0.1*i) for i in 1:Length]
# y=u

#各データのスケーリング
# y=y.-(mean(y)-1.0)
#===================================#


#==================各種パラメータ============#
Length=10000#length(y)#データの総数
Tr=4000#学習用データの長さ
T0=100#過渡時間

Learn_time=T0+Tr#学習に利用する時間
test_time=600#テストに利用するデータの長さ

#各種定数の定義
K=1#入力のデータ次元
L=1#出力のデータ次元
N=600#内部ユニットの個数 Nが小さいと，疎行列の作成ができなくなる(非ゼロの割合で決めているので)
esn=ESN_param(K,N,L)

W_in=esn.W_in

W_res=esn.W_rec

W_back=esn.W_back
#===============================#


#=リザバーの中間層の計算=#
Δt=17.0/Length

Lolenz=RK_N(Lolenz_,Δt,Length,[1,1,1])
y=[]
for i in 1:Length
    push!(y,Lolenz[i][1])
end

rn=zeros(N)#ある時刻でのリザバーベクトル

γ=10.0#大きくすれば学習誤差も下がる
σ=10.0#分岐点

println("γ=",γ,",σ=",σ)

X_ESN=ones(Length,N)#中間層用の格納行列

for n in 1:Length
    X_ESN[n:n,:]=rn#時刻nにおける方程式の解
    yn=y[n]
    # if n==1
    #     local yn=0.0
    # else
    #     local yn=y[n-1]
    # end

    local tn=n*Δt

    local K1=γ*Δt*(-rn+tanh.(W_res*rn+σ*W_back*yn))

    local K2=γ*Δt*((-rn+K1/2)+tanh.((W_res*(rn+K1/2))+σ*(W_back*yn)))

    local K3=γ*Δt*(-(rn+K2/2)+tanh.((W_res*(rn+K2/2))+σ*(W_back*yn)))

    local K4=γ*Δt*(-(rn+K3)+tanh.((W_res*(rn+K3))+σ*(W_back*yn)))

    global rn+=(K1+2.0*K2+2.0*K3+K4)/6.0

end


Learn_M=Matrix(hcat(ones(Length),X_ESN))#バイアス項の追加

λ=10^-6
Wout=Learning(Learn_M,y,λ,T0,Learn_time)


err=NMSE(y[T0:Learn_time,:],Learn_M[T0:Learn_time,:]*Wout)
println("err=",err)
# plot(Learn_M[T0:Learn_time,:]*Wout,xlim=(T0+100,T0+200))
# plot!(y[T0:Learn_time,:])

#============================予測値の計算==================================#
r_new=X_ESN[Learn_time:Learn_time,:]'#学習の最終時刻のデータをベクトル化

for n in Learn_time+1:Length

    # if n==Learn_time+1
    #     local yn=y[n-1]
    # else
    local yn=hcat(1.0,r_new')*Wout
    # end

    local K1=γ*Δt*(-r_new+tanh.(W_res*r_new+W_back*yn))

    local K2=γ*Δt*(-(r_new+K1/2)+tanh.(W_res*(r_new+K1/2)+σ*W_back*yn))

    local K3=γ*Δt*(-(r_new+K2/2)+tanh.(W_res*(r_new+K2/2)+σ*W_back*yn))

    local K4=γ*Δt*(-(r_new+K3)+tanh.(W_res*(r_new+K3)+σ*W_back*yn))

    global r_new+=(K1+2.0*K2+2.0*K3+K4)/6.0

    X_ESN[n:n,:]=r_new#時刻nにおける方程式の解
end

Pred=Matrix(hcat(ones(Length),X_ESN))
# Pred=X_ESN

# err_NMSE=NMSE(y[Learn_time+1:Learn_time+1+test_time,:],Pred[Learn_time+1:Learn_time+1+test_time,:]*Wout)
# println("Err=",err_NMSE)


# T=Learn_time
# plot(Pred*Wout,label=("ESN"),xlims=(T-100,T+100),ylims=(minimum(y),maximum(y)))
# plot!(y)
# savefig("macky_glass.png")
