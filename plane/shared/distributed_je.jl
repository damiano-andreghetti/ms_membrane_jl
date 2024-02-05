@everywhere begin
	println("now load libraries")
	using Random, JLD
	include("MCLangevin_je.jl")
	println("libraries loaded, now initialize variables")
	Random.seed!(22)
	L=100
	l=zeros(L,L)
	part=zeros(Bool,L,L)
	k0=10.0
	ld=0.28
	sigma=1.0
	ratio=0.0#10^-5 #phi/kD
	Nav=10000
	gamma=50.0
	N=1000000 #should take apporx 10h
	Ne=25
	Ninf=10000
	W=0.0
	println("variables initialized")
	function gen_circular_clust(cx,cy,radius,L)
		sites=[]
		for x in 1:L, y in 1:L
			if sqrt((cx-x)^2+(cy-y)^2)<radius
				push!(sites, [x,y])
			end
		end
		return sites
	end
	for el in gen_circular_clust(50,50,5,L)
		part[el[1],el[2]]=true
	end
	function launcher(par)
		k=par[1]
		r=par[2]
		if r!=100
			part[50,Int(55+r)]=true
		end
		fold="K="*string(round(k, digits=3))*",r="*string(round(r, digits=2))
		println("launching simulation for k=", k," and r=", r)
		output=run_sim(l,part,ld,N,Nav,Ninf,k,k0,W,sigma,fold)
	end
end
#logarithmically ranged values
function logrange(x1, x2, n)
	return [10^y for y in range(log10(x1), log10(x2), length=n)]
end

Ks=round.(logrange(10.0,10000.0,7),digits=3)
Rs=[1,2,3,10]#[,20,30,35,100]
parameters=[]
for i in 1:length(Ks)
	for j in 1:length(Rs)
		push!(parameters,[Ks[i],Rs[j]])
	end
end
@time pmap(launcher,parameters)
