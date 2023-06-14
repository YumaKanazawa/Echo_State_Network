using DataFrames,CSV
# using Random
# using Plots


Re=40#レイノルズ数
#=======中間層と入力に関するデータ列の作成=======#
#動的円柱
# z=Matrix(CSV.read("inputdata/input_move"*string(Re)*"80sin.csv",DataFrame,header=false))#入力がsin波
# X=Matrix(CSV.read("reservoir_data/fluid"*string(Re)*"80sin_move.csv",DataFrame,header=false))#データを取り出して行列に格納する.入力がsin波

# z=Matrix(CSV.read("inputdata/input_move"*string(Re)*"80rand.csv",DataFrame,header=false))#入力が乱数
# X=Matrix(CSV.read("reservoir_data/fluid"*string(Re)*"80rand_move.csv",DataFrame,header=false))#データを取り出して行列に格納する.入力が乱数

#ディリクレ
# z=Matrix(CSV.read("inputdata/input"*string(Re)*"_rand.csv",DataFrame,header=false))#乱数入力
# X=Matrix(CSV.read("reservoir_data/fluid"*string(Re)*"_rand_ans.csv",DataFrame,header=false))#乱数入力

# z=Matrix(CSV.read("inputdata/input"*string(Re)*"_sin.csv",DataFrame,header=false))#sin波入力0.25x
# X=Matrix(CSV.read("reservoir_data/fluid"*string(Re)*"_ans_sin.csv",DataFrame,header=false))#sin波入力

#動的円柱連続データ(マッキーグラス)
# z=Matrix(CSV.read("inputdata/input_macky.csv",DataFrame,header=false))
# X=Matrix(CSV.read("reservoir_data/fluid40_macky.csv",DataFrame,header=false))

#動的円柱連続データ(sin波)
# z=Matrix(CSV.read("inputdata/input_cont_sin_m.csv",DataFrame,header=false))
# X=Matrix(CSV.read("reservoir_data/fluid40_cont_sin_m.csv",DataFrame,header=false))

#動的ディリクレ連続データ
# z=Matrix(CSV.read("inputdata/input_cont_sin_d.csv",DataFrame,header=false))
# X=Matrix(CSV.read("reservoir_data/fluid40_cont_sin_d.csv",DataFrame,header=false))
#==========================================#

#NARMAモデルの構成
# function NARMA(n)
#     y=[0.0 for i in 1:n]#nケ分の初期値
#     for k in n:length(z)-1
#         if n==2
#             y_new=(1/1.5)*(0.4*y[k]+0.4*y[k]*y[k-1]+0.6*(z[k]^3)+0.1)
#         else
#             Sy=[y[k-j] for j in 0:n-1]
#             #y[k]+y[k-1]+...+y[k-n+1]
#             y_new=(1/1.95)*(0.3*y[k]+0.05*y[k]*sum(Sy)+1.5*(z[k-n+1]*z[k])+0.1)
#         end
#         push!(y,y_new)
        
#     end
#     return y
# end
#ans_y=NARMA(2)[100:end]
#Z=[0.5*sin(k*0.2/10) for k in 1:length(z)]
#scatter(Z,marker=",",markersize=0.1,markeredgewidth=0.01,label=("y_k"),xlabel=("k"),ylabel=("y_k"))


#plot!(x->0.5*sin(10x),xlims=(10,20))
#plot(ans_y,xlim=(100,200),marker="o",label=("y_k"),xlabel=("k"),ylabel=("y_k"),ylim=(minimum(ans_y),maximum(ans_y)))