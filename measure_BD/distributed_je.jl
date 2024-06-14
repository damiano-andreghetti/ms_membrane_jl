@everywhere begin
	# println("now load libraries")
	using Random, JLD
	include("MCLxy_no_diff.jl")
	println("libraries loaded, now initialize variables")
	Random.seed!(2)
	L=100
	l=zeros(L,L)
	part=zeros(Bool,L,L)
	k0=10.0
	ld=0.25
	ratio=0.0#10^-5 #phi/kD
	Nav=10000
	gamma=500.0
	N=6000000 #should take apporx 20h
	Ne=25
	Ninf=1000
	W=0.0
	println("variables initialized")
	part[50,50]=true
	function launcher(par)
		k=par
		fold="second_run_K="*string(round(k, digits=3))
		println("launching simulation for k=",k)
		output=run_sim(l,part,ld,N,Nav,Ninf,k,k0,W,ratio, Ne, gamma,fold)
	end
end
#logarithmically ranged values
function logrange(x1, x2, n)
	return [10^y for y in range(log10(x1), log10(x2), length=n)]
end

Ks=[10.0,30.0,100.0,300.0]
parameters=[]
for i in 1:length(Ks)
	push!(parameters,Ks[i])
end
@time pmap(launcher,parameters)
