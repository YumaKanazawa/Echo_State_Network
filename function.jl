using LinearAlgebra,SparseArrays
using Random,Distributions,Statistics

Random.seed!(0)#シードの固定
ρ::Float64=0.9
#=============================関数領域=====================================#

function norm(x::Array{Float64,1},p::Float64)
    ret=0.0
    for i in 1:length(x)
        ret+=(x[i])^p
    end
    return (ret/length(x))^(1/p)
end

# NARMAモデルの構成
function NARMA(n,z)
    Y=zeros(Length,L)

    for l in 1:L
        z_l=z[:,l:l]

        z_l=0.4.*(z_l.-minimum(z_l))./(maximum(z_l)-minimum(z_l))

        y=[0.0 for i in 1:n]#nケ分の初期値
        
        for k in n:length(z_l)-1
            if n==2
                y_new=(0.4*y[k]+0.4*y[k]*y[k-1]+0.6*(z_l[k]^3)+0.1)
            else
                Sy=[y[k-j] for j in 0:n-1]
                #y[k]+y[k-1]+...+y[k-n+1]
                y_new=(0.3*y[k]+0.05*y[k]*sum(Sy)+1.5*(z_l[k-n+1]*z_l[k])+0.1)
            end
            push!(y,y_new)
        end
        Y[:,l:l]=y
    end
    return Y
end

NMSE(y::Array{Float64,1},y_h::Array{Float64,1})=sum((y-y_h).^2)/sum((y).^2)#誤差関数

function activate(x)
    return tanh.(x)
end

function Newton(func)
    h=0.001
    x=0.5#rand()#解の初期値
    df=(func(x+h)-func(x-h))/(2*h)

    while true

        # print(i,",",(func(x)/df),",")
        if abs(func(x)/df)<10^-5
            break
        end
        x_new=x-0.01*func(x)/df
        x=x_new
    end
    return x
end

function binary(f)
    #[a,b]内部での解を見つける
    a=0.0
    b=rand()#初期値
    for i in 1:100
        if f(a)*f(b)<=0.0
            break
        end
        b=b+10^3*i
    end

    ret=0.0
    for n in 1:10^5

        if abs(b-a)<10^-12
            break
        end

        c=(a+b)/2
        if f(c)<0.0
            b=c
        else
            a=c
        end
        ret=c

        # println("a=",a,",b=",b)
        # println("ans=",f(ret))
        # println("f(a)=",f(a),",f(b)=",f(b))
    end
    
    return ret
end

function read_file(file_name)
    f=open(file_name,"r")
    list = readlines(f)
    close(f)
    return list
end


function RK_N(F,Δt,T,yn)#y0がベクトルであると考える．
    Ans=[[0.0 for i in 1:length(yn)] for t in 1:T]
    for n in 1:T
        Ans[n]=yn#時刻nにおける方程式の解
        tn=n*Δt

        K1=Δt*F(tn,yn)

        K2=Δt*F(tn+Δt/2,yn+K1/2)

        K3=Δt*F(tn+Δt/2,yn+K2/2)

        K4=Δt*F(tn+Δt,yn+K3/2)

        yn+=(K1+2*K2+2*K3+K4)/6.0
    
        # push!(Ans,yn)
    end
    return Ans
end

mutable struct ESN_param
    N::Int128
    K::Int128
    L::Int128

    W_in::Matrix{Float64}
    W_in_V::Vector{Float64}
    W_rec::Matrix{Float64}
    W_back::Matrix{Float64}
    W_back_V::Vector{Float64}

    #スパース行列
    function SP_matrix(m,n,p)
        sp=Array(sprand(m,n,p))
        for i in 1:size(sp,1)
            for j in 1:size(sp,2)
                if sp[i,j]!=0.0
                    # sp[i,j]=1
                    sp[i,j]=rand(Uniform(-1,1),1)[1]
                    # sp[i,j]=rand(Bernoulli(0.5),1)[1]
                end
            end
        end
        return sp
    end

    #スパースベクトル
    function SP_vector(n,p)
        sp=Array(sprand(n,p))
        for i in 1:length(sp)
            if sp[i]!=0.0
                # sp[i]=1.0
                sp[i]=rand(Uniform(-1,1),1)[1]
            end
        end
        return sp
    end

    #ESP付与
    function spector(a,W)::Matrix{Float64}
        λ_max=maximum(abs.(eigvals(W)))#最大の固有値を求める
        W_return=(1/λ_max).*W
        return a.*W_return
    end

    #行列の初期化
    function ESN_param(K::Int64,N::Int64,L::Int64)
        self=new()
        self.N=N
        self.K=K
        self.L=L
        self.W_in=SP_matrix(N,K,0.0)
        self.W_rec=spector(ρ,SP_matrix(N,N,0.02))
        self.W_back=ones(N,L)#SP_matrix(N,L,1.0)

        self.W_back_V=rand(Uniform(-1,1),N)#出力が1次元の場合の重みベクトル
        self.W_in_V=sprandn(N,0.0)#入力が1次元の場合の重みベクトル
        return self
    end
end

function Learning(X_Learn::Array{Float64,2},Ans::Array{Float64,1},λ::Float64,T0::Int64,Learn_time::Int64)
    M_learn::Array{Float64,2}=zeros(size(X_Learn,1),size(X_Learn,2))#X_Learn[T0:Learn_time,:]#学習に使うデータを上手に切り出す必要がある．N+1はバイアスを含めた次元
    T_learn::Array{Float64,1}=zeros(size(X_Learn,1))

    #X_learn,Ansをτごとに切り出し
    k=1
    for i in 1:Learn_time*τ
        if i==k*τ#τごとに切り出して学習を行う
            M_learn[k:k,:]=X_Learn[i:i,:]
            T_learn[k]=Ans[i]
            k+=1
        end
    end

    use_learn=M_learn[T0:Learn_time,:]
    use_T=T_learn[T0:Learn_time]

    W_out=inv(use_learn'*use_learn+λ*I)*use_learn'*use_T
    
    return W_out
end


#==========================Akaike Infomation Critism(λ's determine method)====================================#
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


# λ=0.001#正規化パラメータ
# # df(λ)の整数解の計算
# AIC_0=10^4#閾値
# for n in 1:N
#     f(x)=df(x)-n
#     # λ_new=Newton(f)#Newton法によるdfの整数解
#     λ_new=binary(f)#方程式の解
#     AIC_new=AIC(λ_new)
#     # println(n,",",df(λ_new),",",AIC_new)
#     # println("λ=",λ,",λ_new=",λ_new)
#     # println()
#     if AIC_0>AIC_new
#         global λ=λ_new
#         global AIC_0=AIC_new
#     end
# end

#========================================================================================#


p=10.0
r=28.0
b=8.0/3.0

function Lolenz_(t,xn)
    #x...とすると行列として定義される
    x=xn
    ret=[0.0 for i in 1:3]
    ret[1]=-p*x[1]+p*x[2]
    ret[2]=-x[1]*x[3]+r*x[1]-x[2]
    ret[3]=x[1]*x[2]-b*x[3]
    return ret
end

function bane(t,xn...)
    x=xn[1]
    ret=[0.0 for i in 1:2]
    ret[1]=-x[2]
    ret[2]=x[1]
    return ret
end


function δ(i::Int64,j::Int64)
    ret=0.0;
    if i==j
        ret=1.0
    else
        ret=0.0
    end
    return ret
end