#内部状態を力学型に変更する
using LinearAlgebra 
using Plots

# NARMAモデルの構成
function NARMA(n,z)
    y=[0.0 for i in 1:n]#nケ分の初期値
    for k in n:length(z)-1
        if n==2
            y_new=(1/1.5)*(0.4*y[k]+0.4*y[k]*y[k-1]+0.6*(z[k]^3)+0.1)
        else
            Sy=[y[k-j] for j in 0:n-1]
            #y[k]+y[k-1]+...+y[k-n+1]
            y_new=(1/1.95)*(0.3*y[k]+0.05*y[k]*sum(Sy)+1.5*(z[k-n+1]*z[k])+0.1)
        end
        push!(y,y_new)
        
    end
    return y
end

Length=2100#とったデータの総数
Tr=2000#学習用データの総数
T0=100#過渡時間
Learn_time=T0+Tr#学習に利用する時間
test_time=Length-Learn_time#テストに利用する長さ

#=============正解データ==============#
max=1.01
min=0.99
Z=[rand([min,max]+rand() for i in Length)]
n=2
Y=NARMA(n,Z)

#===================================#

plot_error=[]

for l in 0:length(Y_mat)-1
#l=50#周波数wと番号lはl=100(w-0.1)という関係
    Y=Y_mat[l+1]
    #========教師収集行列の作成,学習========#
    #学習に利用するもの．100~2100まで

    M=X[T0:Learn_time, :]#T0以降の行を取り出す
    T=Y[T0:Learn_time]#T0以降のデータを抽出する

    #最適な重みベクトルの設定
    λ=10^-5
    W_out=inv(M'*M+λ*I)*M'*T
    ans=X*W_out
    #===================================#

    #========テスト結果の誤差計算========#
    X_test=X[Learn_time+1:end,:]
    T_test=Y[Learn_time+1:end,:]
    predict_error=sum((X_test*W_out-T_test).^2)#
    #===================================#
    
    #==================結果の表示=================#
    #if predict_error<1 && sum((M*W_out-T).^2)<1
        println("周波数=",start+step*l)
        println("予測誤差=",predict_error)
        println("学習誤差=",sum((M*W_out-T).^2))
        println()
    #end
    #グラフの描画
    # Range=1600
    # plot(Y,xlabel=("k"),ylabel=("y^,y"),label=("y"),title=("Re="*string(Re)*"(1/2)sin(x/10)"),xlim=(Range,Range+500),ylim=(minimum(T)-0.01,maximum(T)+0.01))#予測データ
    # plot!(ans,label=("y^"))#正解データ
    #============================================#
end
#annotate!(Length*0.7,minimum(T)-0.005,Plots.text("error≈"*string(error)))#注釈の挿入 