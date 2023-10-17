using Random,Plots,JLD, StatsBase
#using ProfileView,StatsBase, LsqFit, LaTeXStrings
gr()
include("MCLan.jl")
Random.seed!(12)
L=100
global l=zeros(L,L)
global part=zeros(Bool,L,L)
Npart=1#Int(floor(L*L*0.3))
while sum(part)<Npart
	part[rand(1:L),rand(1:L)]=true
end
#part[25:35,25:35].=true
k0=10.0
ld=0.215
k=10.0
N=5000
kI=0.0
Ne=10000
W=0.0
gamma=100.0
Ninf=100
fold="test"
@time output=run_sim(l,part,ld,1,Ninf,k,k0,W,gamma,kI,Ne,fold)
@time output=run_sim(l,part,ld,N,Ninf,k,k0,W,gamma,kI,Ne,fold)
display(plot(output[:,7]))
println(var(output[:,7]))
