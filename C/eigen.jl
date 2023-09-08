using LinearAlgebra
using CSV,DataFrames

i=100
M=Matrix(CSV.read("Matrix.csv",DataFrame,header=false))
A=M[1:i,1:i]

maximum(abs.(eigvals(A)))