using Random, JLD
#TODO
#add compression of output files?
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
	return r::Vector{Vector{Int}}=[[PBC(x+1,L),y], [PBC(x-1,L),y], [x,PBC(y+1,L)], [x,PBC(y-1,L)]]
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
function ev_en(l::Matrix{Float64},part::Matrix{Bool},occupied_neighbours::Matrix{Int},curv::Array{Vector{Float64}},x::Int,y::Int,k::Float64,k0::Float64,W::Float64)
	helfrich_contrib::Float64=0.0
	cs=curv[x,y]
	if part[x,y]
		helfrich_contrib=k*(cs[1]^2+cs[2]^2)-W*occupied_neighbours[x,y]
	else
		helfrich_contrib=k0*((cs[1]+cs[2])^2)
	end	
	#r::Float64=(helfrich_contrib)*0.5
	return helfrich_contrib*0.5
end
function ev_curv(l::Matrix{Float64},x::Int,y::Int,L::Int)
	#we assume in the following that lattice spacing a=1
	x1::Int=PBC(x+1,L)
	x1_::Int=PBC(x-1,L)
	y1::Int=PBC(y-1,L)
	y1_::Int=PBC(y-1,L)
	c::Vector{Float64}=[l[x1,y]+l[PBC(x1_,L),y]-2.0*l[x,y],l[x,PBC(y1,L)]+l[x,PBC(y1_,L)]-2.0*l[x,y]]
	return c #c is a vector composed of [cx,cy]
end
function updateRkl(k::Vector{Int},dir::Int,occupied_neighbours::Matrix{Int},curv::Array{Vector{Float64}},kap::Float64,k0::Float64,W::Float64,L::Int)
	#this function should return f_kl=H(n_k=1,n_l=0)-H(n_k=0,n_l=1)
	#also here assume a=1
	f::Float64=0.0
	if dir==1 #x direction
		cplus=curv[PBC(k[1]+1,L),k[2]]
		cmin=curv[PBC(k[1]-1,L),k[2]]
		f=-W*(occupied_neighbours[PBC(k[1]-1,L),k[2]]-occupied_neighbours[PBC(k[1]+1,L),k[2]])+kap*0.5*(cmin[1]^2+cmin[2]^2-cplus[1]^2-cplus[2]^2)+k0*0.5((cplus[1]+cplus[2])^2-(cmin[1]+cmin[2])^2)
	else #y direction
		cplus=curv[k[1],PBC(k[2]+1,L)]
		cmin=curv[k[1],PBC(k[2]-1,L)]
		f=-W*(occupied_neighbours[k[1],PBC(k[2]-1,L)]-occupied_neighbours[k[1],PBC(k[2]+1,L)])+kap*0.5*(cmin[1]^2+cmin[2]^2-cplus[1]^2-cplus[2]^2)+k0*0.5((cplus[1]+cplus[2])^2-(cmin[1]+cmin[2])^2)
	end
	return f #/(kap-k0)
end
function run_sim(l::Matrix{Float64}, part::Matrix{Bool}, ld::Float64, N::Int,Ninf::Int,kap::Float64,k0::Float64,W::Float64,gamma::Float64,kI::Float64,Ne::Int,folder::String)
	L=Int(sqrt(length(part))) #this works for square lattices only
	measures=zeros(N,7)
	cnt_f=0
	#setup matrix of neighbours
	neig=fill(Int[0,0],(L,L,4))
	Rkl=zeros(Float64,L,L,2)
	for x in 1:L, y in 1:L
		nn=neigb(x,y,L)
		for i in 1:4
			neig[x,y,i]=[nn[i][1],nn[i][2]]
		end
	end
	#setup matrix of number of occupied neighbours
	occ_neig=zeros(Int,L,L)
	for x in 1:L, y in 1:L
		occ_n::Int=0
		@views for el in neig[x,y,:]
			if part[el[1],el[2]]
				occ_n+=1
			end
		end
		occ_neig[x,y]=occ_n
	end
	#setup matrix of curvatures
	curvatures=fill(Float64[0.0,0.0], (L,L))
	for x in 1:L, y in 1:L
		curvatures[x,y]=ev_curv(l,x,y,L)
	end
	#setup matrix of energies per siite
	Ene=zeros(L,L)
	for x in 1:L,y in 1:L
		Ene[x,y]=ev_en(l,part,occ_neig,curvatures,x,y,kap,k0,W)
	end
	for sweep_MC in 1:N
		ord=shuffle(collect(1:L*L))
		acc=0
		free_part_count=0
		acc_free_jump=0
		for st_memb in ord
			#displaccement move
			x,y=num_to_ind(st_memb,L)
			l0=deepcopy(l[x,y])
			Ein::Float64=0.0
			@views for el in neig[x,y,1:4]
				Ein+=Ene[el[1],el[2]]
			end
			Ein+=Ene[x,y]
			delta=((rand()-0.5)*2*ld)	
			l[x,y]+=delta
			@views for el in neig[x,y,1:4]
				curvatures[el[1],el[2]]=ev_curv(l,el[1],el[2],L)
			end
			curvatures[x,y]=ev_curv(l,x,y,L)	
			Efin::Float64=0.0
			@views for el in neig[x,y,1:4]
				Ene[el[1],el[2]]=ev_en(l,part,occ_neig,curvatures,el[1],el[2],kap,k0,W)
				Efin+=Ene[el[1],el[2]]
			end
			Ene[x,y]=ev_en(l,part,occ_neig,curvatures,x,y,kap,k0,W)
			Efin+=Ene[x,y]
			if rand() < exp(-(Efin-Ein))
				#accept
				acc+=1
			else
				#reject and restore config, first l0, then curvatures then energies
				l[x,y]=l0
				@views for el in neig[x,y,1:4]
					curvatures[el[1],el[2]]=ev_curv(l,el[1],el[2],L)
				end
				curvatures[x,y]=ev_curv(l,x,y,L)
				@views for el in neig[x,y,1:4]
					Ene[el[1],el[2]]=ev_en(l,part,occ_neig,curvatures,el[1],el[2],kap,k0,W)
				end
				Ene[x,y]=ev_en(l,part,occ_neig,curvatures,x,y,kap,k0,W)
			end
		end
		#now sweep for particles update all Fkl
		ord=shuffle(collect(1:L*L))
		for st_part in ord
			x,y=num_to_ind(st_part,L)
			if part[x,y]
				for direction in 1:2 #directions x and y
					f=updateRkl([x,y],direction,occ_neig,curvatures,kap,k0,W,L)
					measures[sweep_MC,7]+=f
					cnt_f+=1
					noise=sqrt(2*gamma)*randn()
					#println("f=",f," noise=", noise)
					Rkl[x,y,direction]+=(f+noise)/gamma
				end
			end
		end
		#now check for jump condition
		ord=shuffle(collect(1:L*L))
		for st_part in ord
			x,y=num_to_ind(st_part,L)
			if part[x,y]
				free_part_count+=Int(occ_neig[x,y]==0)
				ab=abs.(Rkl[x,y,1:2])
				if maximum(ab)>0.5
					arr_neig=argmax(ab)
					#idx=0
					if arr_neig==1 && Rkl[x,y,1]>0
						arrival=[PBC(x+1,L),y]
					elseif arr_neig==1 && Rkl[x,y,1]<0
						arrival=[PBC(x-1,L),y]
					elseif arr_neig==2 && Rkl[x,y,2]>0
						arrival=[x,PBC(y+1,L)]
					else
						arrival=[x,PBC(y-1,L)]
					end
					if !part[arrival[1],arrival[2]] #ensure particle conservation
						#bin[idx]+=1
						#jump particles
						part[x,y]=false
						part[arrival[1],arrival[2]]=true
						#update occ_neig
						@views for el in neig[x,y,1:4]
							occ_neig[el[1],el[2]]-=1
						end
						@views for el in neig[arrival[1],arrival[2],1:4]
							occ_neig[el[1],el[2]]+=1
						end
						#update energy
						@views for el in neig[x,y,1:4]
							Ene[el[1],el[2]]=ev_en(l,part,occ_neig,curvatures,el[1],el[2],kap,k0,W)
						end
						@views for el in neig[arrival[1],arrival[2],1:4]
							Ene[el[1],el[2]]=ev_en(l,part,occ_neig,curvatures,el[1],el[2],kap,k0,W)
						end	
						#reset Fkl
						Rkl[x,y,:].=0.0
						Rkl[arrival[1],arrival[2],:].=0.0
						if occ_neig[x,y]==1 && occ_neig[arrival[1],arrival[2]]==0 #if particle is "free"
							acc_free_jump+=1
						end
					else
						#if there's a particle at neighbouring site just stay along the border, without crossing it
						for i in 1:2
							if abs(Rkl[x,y,i])>0.5
								Rkl[x,y,i]=sign(Rkl[x,y,i])*0.5
							end
						end
					end
				end
			end
		end
		#insertion
		ord=shuffle(collect(1:L*L))
		for st_ins in ord
			x,y=num_to_ind(st_ins,L)
			if !part[x,y]
				if rand() < kI
					part[x,y]=true
				end				
			end
		end
		#extraction
		clusters=find_clusters(neig, part, L)
		for cl in clusters
			if length(cl)>=Ne
				#remove particles inside cluster
				for el in cl
					part[el[1],el[2]]=false
				end
			end
		end
		measures[sweep_MC,1]=sum(Ene)
		measures[sweep_MC,2]=acc/(L*L)
		measures[sweep_MC,3]=acc_free_jump/free_part_count
		acc=0
		acc_free_jump=0
		f_nn=0
		interface_size=0
		for x in 1:L, y in 1:L
			@views for nn in neig[x,y,:]
				f_nn+=part[x,y]*part[nn[1],nn[2]]
				interface_size+=part[x,y]*(1-part[nn[1],nn[2]])
			end
		end
		measures[sweep_MC,4]=f_nn/(2*L*L*2)
		measures[sweep_MC,5]=interface_size
		measures[sweep_MC,6]=sum(part)/(L*L)
		measures[sweep_MC,7]/=cnt_f
		cnt_f=0
		if sweep_MC%Ninf==0
			println("sweep ", sweep_MC, " out of ", N)
			println("saving configuration")
			save(folder*"/config_"*string(sweep_MC)*".jld", "particles",part, "membrane", l)
			println("saving measures")
			save(folder*"/measures.jld", "energy",measures[1:sweep_MC,1],"acc_MC",measures[1:sweep_MC,2],"acc_free",measures[1:sweep_MC,3],"f_nn",measures[1:sweep_MC,4],"interface_size",measures[1:sweep_MC,5],"density",measures[1:sweep_MC,6])
			#the following is useful for debugging but MUST BE COMMENTED for cluster usage
			display(heatmap(part, aspect_ratio=:equal,grid=false, c=cgrad(:grays, rev=true), showaxis=false))
		end
	end
	return measures
end

