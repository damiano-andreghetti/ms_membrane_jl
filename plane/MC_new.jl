using Random, JLD
function PBC(x::Int,L::Int)
	rx::Int=0
	if x>L
		rx=x-L
	elseif x<1
		rx=L-x
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
function ev_en(l::Matrix{Float64},part::Matrix{Bool},neig::Array{Vector{Int}},x::Int,y::Int,L::Int,k::Float64,k0::Float64,W::Float64)
	occ_n::Int=0
	helfrich_contrib::Float64=0.0
	cx::Float64=l[PBC(x+1,L),y]+l[PBC(x-1,L),y]-2*l[x,y]
	cy::Float64=l[x,PBC(y+1,L)]+l[x,PBC(y-1,L)]-2*l[x,y]
	if part[x,y]
		helfrich_contrib=k*(cx^2+cy^2)
		@views for el in neig[x,y,:]
			if part[el[1],el[2]]
				occ_n+=1
			end
		end
	else
		helfrich_contrib=k0*((cx+cy)^2)
	end	
	r::Float64=0.5*(helfrich_contrib-W*occ_n)
	return r
end
function membrane_lg(l::Matrix{Float64},part::Matrix{Bool},ld::Float64,N::Int,Nmemb::Int,Nkaw::Int,Ninf::Int,k::Float64,k0::Float64,W::Float64,kI::Float64,Ne::Int,folder::String)
	mis=zeros(N)
	mis2=zeros(N)
	mis3=zeros(N)
	mis_acc_memb=zeros(N)
	mis_acc_kaw=zeros(N)
	mis_eff_ins=zeros(N)
	L=Int(sqrt(length(part))) #this works for square lattices only
	neig=fill(Int[0,0],(L,L,4))
	for x in 1:L, y in 1:L
		nn=neigb(x,y,L)
		for i in 1:4
			neig[x,y,i]=[nn[i][1],nn[i][2]]
		end
	end
	#ld=0.68
	Ene=zeros(L,L)
	for x in 1:L,y in 1:L
		Ene[x,y]=ev_en(l,part,neig,x,y,L,k,k0,W)
	end
	for st in 1:N
		acc_m::Int=0
		acc_k::Int=0
		cont::Int=0
		for sw in 1:Nmemb
			ord=shuffle(collect(1:L*L))
			for st_memb in ord
				#displaccement move
				x,y=num_to_ind(st_memb,L)
				Ein::Float64=0.0
				@views for el in neig[x,y,1:4]
					Ein+=Ene[el[1],el[2]]
				end
				Ein+=Ene[x,y]
				l0=deepcopy(l[x,y])
				delta=((rand()-0.5)*2*ld)	
				l[x,y]+=delta
				Efin::Float64=0.0
				@views for el in neig[x,y,1:4]
					Ene[el[1],el[2]]=ev_en(l,part,neig,el[1],el[2],L,k,k0,W)
					Efin+=Ene[el[1],el[2]]
				end
				Ene[x,y]=ev_en(l,part,neig,x,y,L,k,k0,W)
				Efin+=Ene[x,y]
				p=exp(-(Efin-Ein))
				if rand() < p
					#accept
					acc_m+=1
				else
					l[x,y]=l0
					@views for el in neig[x,y,:]
						Ene[el[1],el[2]]=ev_en(l,part,neig,el[1],el[2],L,k,k0,W)
					end
					Ene[x,y]=ev_en(l,part,neig,x,y,L,k,k0,W)
				end
			end
		end
		mis_acc_memb[st]=acc_m/(L*L*Nmemb)
		#particle diffusion move (Kawasaki)
		for sw in 1:Nkaw
			ord=shuffle(collect(1:L*L))
			for st_kaw in ord
				#displaccement move
				x,y=num_to_ind(st_kaw,L)		
				if part[x,y]
					avail=zeros(Bool, 4)
					for j in 1:4
						nex=neig[x,y,j]
						if !part[nex[1],nex[2]]
							avail[j]=true
						end
					end
					if sum(avail)>0
						free=false
						nex=neig[x,y,shuffle(findall(avail))[1]]
						x2=nex[1]
						y2=nex[2]
						avail2=zeros(Bool, 4)
						for j in 1:4
							nex=neig[x2,y2,j]
							if !part[nex[1],nex[2]]
								avail2[j]=true
							end
						end
						if sum(avail)==4 && sum(avail2)==3 #3 because before move arrival site has one occupied neighbour (departure)
							cont+=1
							free=true
						end
						Ein::Float64=0.0
						@views for el in neig[x,y,:]
							Ein+=Ene[el[1],el[2]]
						end
						@views for el in neig[x2,y2,1:4]
							Ein+=Ene[el[1],el[2]]
						end
						#switch them
						part[x,y]=Bool(1-part[x,y])
						part[x2,y2]=Bool(1-part[x2,y2])
						Efin::Float64=0.0
						@views for el in neig[x,y,:]
							Ene[el[1],el[2]]=ev_en(l,part,neig,el[1],el[2],L,k,k0,W)
							Efin+=Ene[el[1],el[2]]
						end
						@views for el in neig[x2,y2,1:4]
							Ene[el[1],el[2]]=ev_en(l,part,neig,el[1],el[2],L,k,k0,W)
							Efin+=Ene[el[1],el[2]]
						end
						p=exp(-(Efin-Ein))
						if rand() < p
							#accept
							if free
								acc_k+=1
							end
						else
							part[x,y]=Bool(1-part[x,y])
							part[x2,y2]=Bool(1-part[x2,y2])
							@views for el in neig[x,y,:]
								Ene[el[1],el[2]]=ev_en(l,part,neig,el[1],el[2],L,k,k0,W)
							end
							@views for el in neig[x2,y2,:]
								Ene[el[1],el[2]]=ev_en(l,part,neig,el[1],el[2],L,k,k0,W)
							end
						end
					end
				end			
			end
		end
		mis_acc_kaw[st]=acc_k/cont
		mis2[st]=sum(Ene)
		#try insertion
		cont=0
		ord=shuffle(collect(1:L*L))
		for site in ord
			x,y=num_to_ind(site,L)
			if !part[x,y]
				if rand() < kI
					part[x,y]=true
					cont+=1
				end
			end
		end
		mis_eff_ins[st]=cont/(L*L)
		#extraction
		cl=find_clusters(neig,part,L)
		for clust in cl
			if length(clust)>=Ne
				#extract
				for el in clust
					part[el[1],el[2]]=false
				end
			end
		end		
		mis[st]=sum(part)/(L*L)
		corr=0
		for x in 1:L, y in 1:L
			@views for nn in neig[x,y,:]
				corr+=part[x,y]*part[nn[1],nn[2]]
			end
		end
		mis3[st]=corr/(2*L*L*2)
		if st%Ninf==0
			println("completed sweep ", st)
			println("saving configuration")
			save(folder*"config_"*string(st)*".jld", "particles",part, "membrane", l)
			println("saving measures")
			save(folder*"measures.jld", "density", mis[1:st],"energy",mis2[1:st],"correlation_nn",mis3[1:st],"acc_mem",mis_acc_memb[1:st], "acc_kaw",mis_acc_kaw[1:st],"eff_ins",mis_eff_ins[1:st])
			#mv(folder*"measures_new.jld", folder*"measures.jld", force=true) #this is such that if software stops while saving doesnt actually looses data
		end
	end
	return mis, mis2,mis3,mis_acc_memb,mis_acc_kaw,mis_eff_ins #return density, energy, correlation nn, accepteance for membrane moves, acceptance for free particles diffusion moves
end
