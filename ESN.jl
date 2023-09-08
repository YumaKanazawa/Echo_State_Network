include("ESN_setting.jl")
using Plots
#==========================データの学習==========================#

Length::Int64=Learn_time+test_time

# y::Array{Float64,1}=tryparse.(Float64, read_file("Basic_ESN/data/MackeyGlass_t17.txt"))
# y=[sin(0.01*i) for i in 1:Length]


#ネットワークの初期化
X_ESN=zeros(Length,N)#学習に利用する内部状態の行列
#標本抽出(内部状態の作成)
x=zeros(N)'#x(0),横ベクトルで定義する．力学系の初期値として決める

for l in 1:Learn_time#l行目の内部状態を作成する
    local y_l=0.0
    if l==1
        y_l=0.0
    else
        y_l=y[l-1]
    end
    X_ESN[l:l,:]=tanh.(W_res*x'+W_back*y_l)#X[n:n,:]でn行目からn行目まで取り出すことができる
    global x=X_ESN[l:l,:]#x(n+1)=x(n)に更新する
end

Learn=hcat(ones(Length),X_ESN,X_ESN.^2)#バイアス項の追加
λ=10^-6#正規化パラメータ

#==========================学習=================================#
W_out_ESN=Learning(Learn,y,λ,T0,Learn_time)
#==============================================================#

err=NMSE(y[T0:Learn_time],Learn[T0:Learn_time,:]*W_out_ESN)

#============================予測値の計算==================================#
x_pred=X_ESN[Learn_time:Learn_time,:]#学習の最終時刻のデータをベクトル化

for l in Learn_time+1:Length#l行目の内部状態を作成する
    println(l)
    # local u_l=u[l]#u_lは時刻lの入力データベクトル
    if l==Learn_time+1
        local y_l=y[l-1]#各時刻でL次元のベクトル
    else
        local one_past_predict=hcat(1.0,x_pred,x_pred.^2)#１時刻前の予測結果
        local y_l=(one_past_predict*W_out_ESN)[1]#Closed Loop
        # local y_l=y[l-1]#OpenLoop
    end
    X_ESN[l:l,:]=tanh.(W_res*x_pred'+W_back*y_l)#X[n:n,:]でn行目からn行目まで取り出すことができる
    global x_pred=X_ESN[l:l,:]#x(n+1)=x(n)に更新する
end

println("Train err=",err)

Pred=Matrix(hcat(ones(Length),X_ESN,X_ESN.^2))
# Pred=X_ESN

err_NMSE=NMSE(y[Learn_time+1:Learn_time+test_time],(Pred[Learn_time+1:Learn_time+test_time,:]*W_out_ESN))
println("Err=",err_NMSE)
T=Learn_time
plot(Pred*W_out_ESN,label=("ESN"),xlims=(T,Learn_time+test_time),ylims=(minimum(y),maximum(y)))
plot!(y)
# savefig("macky_glass.png")
