using Random, JLD
#todo
#add insertion with kI such that phi/kD is fixed
#add extraction

function PBC(x::Int,L::Int)
	rx::Int=0
	if x>L
		rx=x-L
	elseif x<1
		rx=L+x
	else
		rx=x
	end
	return rx
end
function num_to_ind(num::Int, L::Int)
	return ((num-1)%L)+1, div(num-1, L)+1
end
function neigb(x::Int,y::Int,L::Int)
	return r::Vector{Vector{Int}}=[[PBC(x+1,L),y], [PBC(x-1,L),y], [x,PBC(y+1,L)], [x,PBC(y-1,L)], [PBC(x-1,L),PBC(y-1,L)],[PBC(x-1,L),PBC(y+1,L)],[PBC(x+1,L),PBC(y+1,L)],[PBC(x+1,L),PBC(y-1,L)]]
end
#function neigb_2(x::Int,y::Int,L::Int)
#	return r::Vector{Vector{Int}}=[[PBC(x+1,L),PBC(y+1,L)], [PBC(x+2,L),y], [PBC(x+1,L),PBC(y-1,L)],[x,PBC(y-2,L)], [PBC(x-1,L),PBC(y-1,L)],[PBC(x-2,L),y],[PBC(x-1,L),PBC(y+1,L)],[x,PBC(y+2,L)]]
#end
function find_clusters(neig::Array{Vector{Int},3},part::Matrix{Bool},L::Int)
	function add_ignore(x,y, ign)
		clust=ign
		@views for j in neig[x,y,1:4]
			if part[j[1],j[2]] && !(j in clust)
				push!(clust, j)
				clust=add_ignore(j[1],j[2], clust)
			end
		end
		return clust
	end
	visited=zeros(Bool,L,L)
	clusters=[]
	for x in 1:L, y in 1:L
		if part[x,y] && !visited[x,y]		
			cc=add_ignore(x,y,[[x,y]])
			push!(clusters, cc)
			for el in cc
				visited[el[1],el[2]]=true
			end
		end
	end
	return clusters
end
function ev_en(l::Matrix{Float64},part::Matrix{Bool},occupied_neighbours::Matrix{Int},derivatives::Array{Vector{Float64}},x::Int,y::Int,k::Float64,k0::Float64,W::Float64,sigma::Float64)
	ds=derivatives[x,y]
	dA::Float64=sqrt(1+ds[3]^2+ds[4]^2)
	M::Float64=(ds[1]*(1+ds[4]^2)+ds[2]*(1+ds[3]^2)-2*ds[5]*ds[3]*ds[4])/(dA^3) #this is actually 2*M = (c1+c2)
	G::Float64=(ds[1]*ds[2]-ds[5]^2)/(dA^4)
	energy::Float64=(M)^2 -2*G #(c1+c2)^2 -2*c1*c2
	if part[x,y]
		energy=k*energy
	else
		energy=k0*energy
	end	
	return energy*0.5*dA+log(dA)+sigma*dA-(W*occupied_neighbours[x,y]*part[x,y]*0.5)
end
function ev_der(l::Matrix{Float64},x::Int,y::Int,L::Int)
	#we assume in the following that lattice spacing a=1, evaluate derivatives, in order h_xx, h_yy, h_x, h_y, h_xy7
	x1::Int=PBC(x+1,L)
	x1_::Int=PBC(x-1,L)
	y1::Int=PBC(y+1,L)
	y1_::Int=PBC(y-1,L)
	h_xx::Float64=l[x1,y]+l[x1_,y]-2.0*l[x,y]
	h_yy::Float64=l[x,y1]+l[x,y1_]-2.0*l[x,y]
	h_x::Float64=(l[x1,y]-l[x1_,y])/2.0
	h_y::Float64=(l[x,y1]-l[x,y1_])/2.0
	h_xy::Float64=(l[x1,y1]-l[x1_,y1]-l[x1,y1_]+l[x1_,y1_])/4.0
	return [h_xx,h_yy,h_x,h_y,h_xy]
end
function ev_der(l::Matrix{Float64},x::Int,y::Int,L::Int,i::Int)
	#we assume in the following that lattice spacing a=1, evaluate derivatives, in order h_xx, h_yy, h_x, h_y, h_xy7
	x1::Int=PBC(x+1,L)
	x1_::Int=PBC(x-1,L)
	y1::Int=PBC(y+1,L)
	y1_::Int=PBC(y-1,L)
	if i==1
		return h_xx::Float64=l[x1,y]+l[x1_,y]-2.0*l[x,y]
	elseif i==2
		return h_yy::Float64=l[x,y1]+l[x,y1_]-2.0*l[x,y]
	elseif i==3
		return h_x::Float64=(l[x1,y]-l[x1_,y])/2.0
	elseif i==4
		return h_y::Float64=(l[x,y1]-l[x,y1_])/2.0
	elseif i==5
		return h_xy::Float64=(l[x1,y1]-l[x1_,y1]-l[x1,y1_]+l[x1_,y1_])/4.0
	end
end

function comp_DH(a::Int,b::Int,c::Int,d::Int,derivatives::Array{Vector{Float64}},o_n::Matrix{Int},k::Float64,k0::Float64,W::Float64)
	f::Float64=0.0
	ds=derivatives[a,b]
	dA::Float64=sqrt(1+ds[3]^2+ds[4]^2)
	M::Float64=(ds[1]*(1+ds[4]^2)+ds[2]*(1+ds[3]^2)-2*ds[5]*ds[3]*ds[4])/(dA^3) #this is actually 2*M = (c1+c2)
	G::Float64=(ds[1]*ds[2]-ds[5]^2)/(dA^4)
	f+=((M^2 -2*G)*(k-k0)*0.5*dA)
	ds=derivatives[c,d]
	dA=sqrt(1+ds[3]^2+ds[4]^2)
	M=(ds[1]*(1+ds[4]^2)+ds[2]*(1+ds[3]^2)-2*ds[5]*ds[3]*ds[4])/(dA^3) #this is actually 2*M = (c1+c2)
	G=(ds[1]*ds[2]-ds[5]^2)/(dA^4)
	f+=((M^2 -2*G)*(k0-k)*0.5*dA)
	f-=W*(o_n[a,b]-o_n[c,d]+1)
	return f
end




function run_sim(l::Matrix{Float64}, part::Matrix{Bool}, ld::Float64, N::Int,Nav::Int,Ninf::Int,kap::Float64,k0::Float64,W::Float64,sigma::Float64,folder::String)
	L=Int(sqrt(length(part))) #this works for square lattices only
	measures=zeros(N,2)
	kds=zeros(N)
	Rkl=zeros(Float64,L,L,4)
	#setup matrix of neighbours
	neig=fill(Int[0,0],(L,L,8))
	#for the membrane part all 8 neighbours are needed, since changing the height at a site involves all these neighbours energies
	#for the particles part just the first 4 are needed, they just go up, down, left, right
	for x in 1:L, y in 1:L
		nn=neigb(x,y,L)
		for i in 1:8
			neig[x,y,i]=[nn[i][1],nn[i][2]]
		end
	end
	#setup matrix of number of occupied neighbours
	occ_neig=zeros(Int,L,L)
	for x in 1:L, y in 1:L
		occ_n::Int=0
		@views for el in neig[x,y,1:4]
			if part[el[1],el[2]]
				occ_n+=1
			end
		end
		occ_neig[x,y]=occ_n
	end
	#setup matrix of curvatures
	derivatives=fill(Float64[0.0,0.0,0.0,0.0,0.0], (L,L))
	for x in 1:L, y in 1:L
		derivatives[x,y]=ev_der(l,x,y,L)
	end
	#setup matrix of energies per siite
	Ene=zeros(L,L)
	for x in 1:L,y in 1:L
		Ene[x,y]=ev_en(l,part,occ_neig,derivatives,x,y,kap,k0,W,sigma)
	end
	for sweep_MC in 1:N
		ord=shuffle(collect(1:L*L))
		acc=0
		free_jump_possible=0
		free_jump_actual=0
		for st_memb in ord
			#displaccement move
			x,y=num_to_ind(st_memb,L)
			l0=deepcopy(l[x,y])
			Ein::Float64=0.0
			@views for el in neig[x,y,1:8]
				Ein+=Ene[el[1],el[2]]
			end
			Ein+=Ene[x,y]
			delta=((rand()-0.5)*2*ld)	
			l[x,y]+=delta
			@views for el in neig[x,y,1:4]
				derivatives[el[1],el[2]]=ev_der(l,el[1],el[2],L)
			end
			@views for el in neig[x,y,5:8]
				derivatives[el[1],el[2]][5]=ev_der(l,el[1],el[2],L,5)
			end
			derivatives[x,y]=ev_der(l,x,y,L)	
			Efin::Float64=0.0
			@views for el in neig[x,y,1:8]
				Ene[el[1],el[2]]=ev_en(l,part,occ_neig,derivatives,el[1],el[2],kap,k0,W,sigma)
				Efin+=Ene[el[1],el[2]]
			end
			Ene[x,y]=ev_en(l,part,occ_neig,derivatives,x,y,kap,k0,W,sigma)
			Efin+=Ene[x,y]
			if rand() < exp(-(Efin-Ein))
				#accept
				acc+=1
			else
				#reject and restore config, first l0, then curvatures then energies
				l[x,y]=l0
				@views for el in neig[x,y,1:4]
					derivatives[el[1],el[2]]=ev_der(l,el[1],el[2],L)
				end
				@views for el in neig[x,y,5:8]
					derivatives[el[1],el[2]][5]=ev_der(l,el[1],el[2],L,5)
				end
				derivatives[x,y]=ev_der(l,x,y,L)	
				@views for el in neig[x,y,1:8]
					Ene[el[1],el[2]]=ev_en(l,part,occ_neig,derivatives,el[1],el[2],kap,k0,W,sigma)
				end
				Ene[x,y]=ev_en(l,part,occ_neig,derivatives,x,y,kap,k0,W,sigma)
			end
		end
		measures[sweep_MC,1]=sum(Ene)
		measures[sweep_MC,2]=acc/(L*L)
		acc=0
		if sweep_MC%Ninf==0
			println("sweep ", sweep_MC, " out of ", N)
			println("saving configuration")
			save(folder*"/config_"*string(sweep_MC)*".jld", "particles",part, "membrane", l, "Rkl", Rkl)
			println("saving measures")
			save(folder*"/measures.jld", "energy",measures[1:sweep_MC,1],"acc_MC",measures[1:sweep_MC,2])
			#the following is useful for debugging but MUST BE COMMENTED for cluster usage
			#display(plot(measures[1:sweep_MC,2]))
			#display(plot(xs,hh))
			#display(plot([heatmap(l, aspect_ratio=:equal),heatmap(part, aspect_ratio=:equal,grid=false, c=cgrad(:grays, rev=true), showaxis=false)]...,layout=(2,1)))
			#display(plot(kds[1:sweep_MC]))
			#display(heatmap(part, aspect_ratio=:equal,grid=false, c=cgrad(:grays, rev=true), showaxis=false))
		end
	end
	return measures
end


