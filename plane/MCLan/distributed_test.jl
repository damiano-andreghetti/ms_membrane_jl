@everywhere begin
	println("now load libraries")
	using Random, JLD
	include("MCLangevin_v3.jl")
	println("libraries loaded, now initialize variables")
	Random.seed!(22)
	L=100
	l=zeros(L,L)
	part=zeros(Bool,L,L)
	k0=10.0
	ld=0.215
	sigma=0.6
	ratio=10^-5 #phi/kD
	Nav=10000
	gamma=50.0
	N=3000000 #should take apporx 20h
	Ne=25
	Ninf=10000
	W=0.0
	println("variables initialized")
	#k=25.0
	println("variables initialized")
	function launcher(par)
		k=par
		fold="K="*string(round(k, digits=3))*",W="*string(round(W, digits=2))
		println("launching simulation for k=", k," and W=", W)
		output=run_sim(l,part,ld,N,Nav,Ninf,k,k0,W,sigma,ratio,Ne,gamma,fold)
	end
end
#logarithmically ranged values
function logrange(x1, x2, n)
	return [10^y for y in range(log10(x1), log10(x2), length=n)]
end

Ks=round.(logrange(10.0,1000.0,9),digits=3)
parameters=[]
for i in 1:length(Ks)
	push!(parameters,Ks[i])
end
@time pmap(launcher,parameters)
