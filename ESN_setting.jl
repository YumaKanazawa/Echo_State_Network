include("function.jl")
using Plots

#=============教師データ==============#
# max=1.01
# min=0.99

# u=rand(Uniform(min,max),DataLength,K)#一様乱数を取得する
# n=2
# y=NARMA(n,u)

# y::Array{Float64,1}=tryparse.(Float64, read_file("../Basic_ESN/data/MackeyGlass_t17.txt"))
# y=u

# u=[sin(0.1*i) for i in 1:DataLength]
# y::Array{Float64,1}=u

#各データのスケーリング
# y=y.-(mean(y)-1.0)
#===================================#


Δt::Float64=0.01#数値計算の刻み幅

#=================各種パラメータ=================#
τ=1#データを取る間隔
Tr::Int64=7000#学習用データの長さ
T0::Int64=100#過渡時間
Learn_time::Int64=T0+Tr#学習に利用する時間
test_time::Int64=10^3#テストに利用するデータの長さ

#=======シミュレーションに利用するデータの長さ=======#
DataLength::Int64=τ*Learn_time+test_time#length(y)#データの総数

# 各種定数の定義
K::Int64=1#入力のデータ次元
L::Int64=1#出力のデータ次元
N::Int64=5000#内部ユニットの個数 Nが小さいと，疎行列の作成ができなくなる(非ゼロの割合で決めているので)
esn=ESN_param(K,N,L)

# W_in::Array{Float64,2}=esn.W_in
W_in::Array{Float64,1}=esn.W_in_V

W_res::Array{Float64,2}=esn.W_rec

W_back::Array{Float64,1}=esn.W_back_V
#===============================#

#=ローレンツアトラクタの計算=#
#=リザバーの中間層の計算=#

Lolenz=RK_N(Lolenz_,Δt,DataLength,[1,1,1])
y::Array{Float64,1}=zeros(DataLength)
for i in 1:DataLength
    # l=rand(1:3)
    y[i]=Lolenz[i][2]
    #y[i]=y[i]/abs(y[i])#大きさを1に正規化
end

γ::Float64=10.0#大きくすれば学習誤差も下がる
σ::Float64=0.012#分岐点

println("γ=",γ,",σ=",σ)
#=================この上までがリザバーの基本的な設定========================#

#ヤコビ行列
function Jacobi_matrix(γ::Float64,σ::Float64,n::Int64,xn::Array{Float64,1})::Array{Float64,2}
    J::Array{Float64,2}=zeros(N,N)
    s=W_res*xn+σ*W_back*y[n]
    for i in 1:N
        for j in 1:N
            J[i,j]=γ*(-δ(i,j)+W_res[i,j]/cosh(s[i]))
        end
    end
    return J
end

# plot(y,xlim=(100,200))