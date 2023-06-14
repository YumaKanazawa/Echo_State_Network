include("function.jl")
using Plots

Length=7100#データの総数
Tr=4000#学習用データの長さ
T0=100#過渡時間

Learn_time=T0+Tr#学習に利用する時間
test_time=2000#テストに利用するデータの長さ

#各種定数の定義
K=1#入力のデータ次元
L=1#出力のデータ次元
N=600#内部ユニットの個数
esn=ESN_param(K,N,L)

W_in=esn.W_in

W_res=esn.W_rec

W_back=esn.W_back

#=============教師データ==============#
# max=1.01
# min=0.99
Random.seed!(0)#シードの固定
# u=rand(Uniform(min,max),Length,K)#一様乱数を取得する
# n=2
# y=NARMA(n,u)

u=tryparse.(Float64, read_file("Basic ESN/data/MackeyGlass_t17.txt"))
y=u

# u=[sin(0.1*i) for i in 1:Length]
# y=u

#各データのスケーリング
# y=y.-(mean(y)-1.0)
#===================================#

#==========================データの学習==========================#

#ネットワークの初期化
X_ESN=zeros(Length,N)#学習に利用する内部状態の行列
#標本抽出(内部状態の作成)
x=zeros(N)'#x(0),横ベクトルで定義する．力学系の初期値として決める

for l in 1:Learn_time#l行目の内部状態を作成する
    local u_l=u[l:l,:]'#u_lは時刻lの入力データベクトル
    if l==1
        local y_l=[0.0 for i in 1:L]#各時刻でL次元のベクトル
    else
        local y_l=y[l-1:l-1,:]'
    end
    X_ESN[l:l,:]=activate(W_in*u_l+W_res*x'+W_back*y_l)#X[n:n,:]でn行目からn行目まで取り出すことができる
    global x=X_ESN[l:l,:]#x(n+1)=x(n)に更新する
end

Learn=Matrix(hcat(ones(Length),X_ESN))#バイアス項の追加
# Learn=X_ESN

# #==========================Akaike Infomation Critism(λ's determine method)====================================#
# function df(λ)
#     AIC_transient=200
#     AIC_length=1200
#     AIC_X=X_ESN[AIC_transient:+AIC_transient+AIC_length,:]#AIC計算のための切り出し
#     # AIC_T=y[AIC_transient:+AIC_transient+AIC_length,:]#AIC計算のための切り出し

#     Lambda=eigvals(AIC_X'*AIC_X)
#     ret=0.0
#     for i in eachindex(Lambda)
#         ret+=Lambda[i]/(abs(λ)+Lambda[i])#正の根を得るために|・|をつける
#     end
#     return ret
# end

# function AIC(λ_func)
#     AIC_transient=200
#     AIC_length=1200
#     AIC_X=X_ESN[AIC_transient:+AIC_transient+AIC_length,:]#AIC計算のための切り出し
#     AIC_T=y[AIC_transient:+AIC_transient+AIC_length,:]#AIC計算のための切り出し

#     Wout=inv(AIC_X'*AIC_X+λ_func*I)*AIC_X'*AIC_T

#     aic=(AIC_length)*log(sum((AIC_X*Wout.-AIC_T).^2))+2*df(λ_func)
#     return aic
# end

# function err(λ_func)
#     AIC_transient=200
#     AIC_length=1200
#     AIC_X=X_ESN[AIC_transient:+AIC_transient+AIC_length,:]#AIC計算のための切り出し
#     AIC_T=y[AIC_transient:+AIC_transient+AIC_length,:]#AIC計算のための切り出し

#     Wout=inv(AIC_X'*AIC_X+λ_func*I)*AIC_X'*AIC_T

#     aic=AIC_length*log(sum((AIC_X*Wout.-AIC_T).^2))
#     return aic
# end


λ=0.001#正規化パラメータ
# # df(λ)の整数解の計算
# AIC_0=10^4#閾値
# for n in 1:N
#     f(x)=df(x)-n
#     # λ_new=Newton(f)#Newton法によるdfの整数解
#     λ_new=binary(f)#方程式の解
#     AIC_new=AIC(λ_new)
#     println(n,",",df(λ_new),",",AIC_new)
#     # println("λ=",λ,",λ_new=",λ_new)
#     # println()
#     if AIC_0>AIC_new
#         global λ=λ_new
#         global AIC_0=AIC_new
#     end
# end
#========================================================================================#


#==========================学習=================================#
W_out_ESN=Learning(Learn,y,λ,T0,Learn_time)
#==============================================================#

# err=NMSE(y[T0:Learn_time,:],Learn[T0:Learn_time,:]*W_out_ESN)
# println("err=",err)
# plot(Learn[T0:Learn_time,:]*W_out_ESN,xlim=(T0+100,T0+200))
# plot!(y[T0:Learn_time,:])

#============================予測値の計算==================================#
x_pred=X_ESN[Learn_time:Learn_time,:]#学習の最終時刻のデータをベクトル化
for l in Learn_time+1:Length#l行目の内部状態を作成する
    local u_l=u[l:l,:]'#u_lは時刻lの入力データベクトル
    if l==Learn_time+1
        local y_l=y[l-1:l-1,:]'#各時刻でL次元のベクトル
    else
        local one_past_predict=hcat(1.0,x_pred)#１時刻前の予測結果
        # local one_past_predict=x_pred#１時刻前の予測結果
        local y_l=one_past_predict*W_out_ESN#この時点では1次元のベクトルとして値を保持
    end
    X_ESN[l:l,:]=activate(W_in*u_l+W_res*x_pred'+W_back*y_l)#X[n:n,:]でn行目からn行目まで取り出すことができる
    global x_pred=X_ESN[l:l,:]#x(n+1)=x(n)に更新する
end

Pred=Matrix(hcat(ones(Length),X_ESN))
# Pred=X_ESN

err_NMSE=NMSE(y[Learn_time+1:Length,:],Pred[Learn_time+1:Length,:]*W_out_ESN)
println("Err=,",err_NMSE)
T=4000
plot(Pred*W_out_ESN,label=("ESN"),xlims=(T,T+200),ylims=(minimum(y),maximum(y)))
plot!(y)
savefig("macky_glass.png")
