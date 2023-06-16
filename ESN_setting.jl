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

#=ローレンツアトラクタの計算=#

Lolenz=RK_N(Lolenz_,0.02,DataLength,[1,1,1])
y::Array{Float64,1}=zeros(DataLength)
for i in 1:DataLength
    y[i]=Lolenz[i][1]
end

rn::Array{Float64,1}=zeros(N)#ある時刻でのリザバーベクトル

γ::Float16=10.0#大きくすれば学習誤差も下がる
σ::Float16=0.01#分岐点

println("γ=",γ,",σ=",σ)

X_ESN::Array{Float64,2}=zeros(DataLength,N)#中間層用の格納行列

#=リザバーの中間層の計算=#
Δt::Float64=0.02#17.0/DataLength

#=================この上までがリザバーの基本的な設定========================#

function δ(i::Int64,j::Int64)
    ret=0.0;
    if i==j
        ret=1.0
    else
        ret=0.0
    end
    return ret
end

#ヤコビ行列
function Jacobi_matrix(γ::Float64,σ::Float64)::Array{Float64,2}
    J::Array{Float64,2}=zeros(N,N)
    s=W_res*rn+σ*W_back*y[DataLength]
    for i in 1:N
        for j in 1:N
            J[i,j]=γ*(δ(i,j)+W_res[i,j]/cosh(s[i]))
        end
    end
    return J
end