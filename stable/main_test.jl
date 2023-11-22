println("now load libraries")
using Random, JLD, Plots
include("MCLangevin.jl")
println("libraries loaded, now initialize variables")
Random.seed!(22)
L=100
l=zeros(L,L)
part=zeros(Bool,L,L)
k0=10.0
ld=0.28
sigma=1.0
ratio=10^-5 #phi/kD
Nav=10000
gamma=50.0
N=500000 #should take apporx 20h
Ne=25
Ninf=1000
#W=0.0
println("variables initialized")
#k=25.0
println("variables initialized")
k=10.0
W=0.0
fold="test_new_version2"
println("launching simulation for k=", k," and W=", W)
output=run_sim(l,part,ld,N,Nav,Ninf,k,k0,W,sigma,ratio,Ne,gamma,fold)
display(plot(output[:,5]))

