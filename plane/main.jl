using Random,Plots, ProfileView,StatsBase, JLD
gr()
include("MC_new.jl")
Random.seed!(1234)
L=100
global l=rand(L,L)
global part=zeros(Bool,L,L)
k0=50.0
ld=0.095
k=2500.0
Nmemb=1 #number of sweeps for membrane, each sweep try to displace in random order each membrane site
Nkaw=0 #number of randomly selected sites, where if there's a particle diffusion is tried
Nrun=2000
Ninf=500 #how often to print step progression and save configuration and measures
W=0.0
kI=0.0#0.00001
Ne=10000
println("stuff loaded")
fold="testk1/"
@time dens,ene,corr,acc_mem,acc_kaw,eff_ins=membrane_lg(l,part,ld,Nrun,Nmemb,Nkaw,Ninf,k,k0,W,kI,Ne,fold)	
plot(acc_mem)