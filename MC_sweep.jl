using Random
#Random.seed!(1234)
include("sphere.jl") #requires TriSurface and CellsGrid definition
I3=MMatrix{3,3,Float64}([1. 0. 0.; 0. 1. 0.; 0. 0. 1.]) #identity 3x3


function ev_energyonset(su::SurfGeoQuant, k::Float64, part::Vector{Int16}, mu::Float64, set::Vector{Int16})
	E::Float64=0
	for n in set
		Av=0
		Nv=MVector{3,Float64}(0.,0.,0.)
		for nf in su.neig_faces[n]
			af=su.face_area[nf]
			Av+=af
			Nv+=(su.face_normals[nf]*af)
			#Nv=sum_vec(Nv,(su.face_normals[nf]*af))
		end
		Nv/=norm(Nv)
		Av/=3
		#Sv=zeros(3,3)
		Sv=MMatrix{3,3,Float64}(0., 0., 0.,0., 0., 0.,0. ,0. ,0.)
		Pv=MMatrix{3,3,Float64}(I3-tensor_prod(Nv,Nv))
		for e in su.neig_edges[n]
			#su.edges[e][1] == n ? j=su.edges[e][2] : j=su.edges[e][1] # j = neighbour of i, with which edge e is in common
			Ne=su.face_normals[su.edge_faces[e][2]]+su.face_normals[su.edge_faces[e][1]]
			Ne=Ne./norm(Ne)
			We=scalar_prod(Nv, Ne)
			Sv=Sv+(We*su.Se[e])
		end
		Sv=Pv*Sv*Pv #notice that transpose(Pv)=Pv and Pv is made of real number => adjoint(Pv)=Pv
		#Householder transformation used to evaluate Sv eigenvalues faster (test against LinearAlgebra.eigvals, is much faster
		#c1,c2=eigenHH(Sv, Nv)
		#println(c1+c2)
		#println(Sv[1,1]+Sv[2,2]+Sv[3,3])
		#cm=(c1+c2)/2
		cm=Sv[1,1]+Sv[2,2]+Sv[3,3]
		"""
		dk=0
		pn=part[n]
		if pn==1 #if a particle is present on vertex, then induce local curvature
			dk=mu
		else
			for j in su.neig[n]
				if part[j]==1
					dk=mu
				end
			end
		end
		t=0.0
		cnt=0
		J=0.0
		for j in su.neig[n]
			cnt+=(pn*part[j])
		end
		E+=((cm^2)*Av*(k+dk)-J*cnt+t*Av)
		"""
		E+=((cm*cm)/Av)
	end
	return E*k
end

#function to evaluate total energy, OLD NEEDS TO BE FIXED
function ev_energy(su::SurfGeoQuant, k::Float64, part::Vector{Int16}, mu::Float64)
	En=ev_energyonset(su,k*0.5,part,mu, [Int16(i) for i in 1:length(su.vertices)])
	return En
end


#evaluate energy difference in case of moving vertex i to new position xn (assumin this is an acceptable position)
function ev_envm(su::SurfGeoQuant, i::Int16, xn::MVector{3,Float64}, k::Float64, part::Vector{Int16}, mu::Float64) 
	#when vertex i is displaced to evaluate energy difference we need to evaluate change in curvatures for i and all his neighbours
	inv_ver=convert(Vector{Int16}, vcat(su.neig[i], [i]))
	Ein=ev_energyonset(su, k, part, mu, inv_ver)
	su2 = update_geoquantvm(su, i, xn)
	Efin=ev_energyonset(su2,k,part,mu, inv_ver)
	return Efin-Ein, su2
end

function ev_enlm(su::SurfGeoQuant, e::Int16, k::Float64, part::Vector{Int16}, mu::Float64, set::Vector{Int16})
	t1,t2=su.edge_faces[e][1],su.edge_faces[e][2]
	Ein=ev_energyonset(su,k,part,mu,set)
	su2, g = update_geoquantlm(su, e, set)
	if g==false
		return 0., su2, g
	end
	Efin=ev_energyonset(su2,k,part,mu,set)
	return Efin-Ein, su2, g
end

#function MC sweep accounting for particles move
function MC_sweepPM(surf::SurfGeoQuant, CELLS::CellsGrid, sigma::Float64, max_edge_len::Float64, site_radius::Float64, k::Float64, part::Vector{Int16}, mu::Float64,Nsvert::Int64, Nslink::Int64, Nspart::Int64; verbose::Bool=false)
	b=1.
	k=k*0.5 #from now on we'll just need k/2 so there's no matter of repeating this operation many times
	a=[0.,0.,0.]
	#mc step for every vertex
	acc=0
	for stepvertex in 1:Nsvert
		vert=Int16(rand(1:surf.Nv))
		x0=deepcopy(surf.vertices[vert])
		#xn=sum_vec(x0,((rand(3).-0.5)*2*sigma))
		delta=MVector{3,Float64}((rand(3).-0.5)*2*sigma)
		xn=x0+delta
		good=true
		#check for max edge distance
		for n in surf.neig[vert]
			if good==false
				break;
			end
			if distance(xn, surf.vertices[n])>max_edge_len
				good=false
			end
			#distance(xn, surf.vertices[n])>max_edge_len && good==true ? good=false : good=true
		end
		if good #check for non overlapping of spheres
			#println("vertex: ", vert, "reasonably around vertices are: ", reasonablyaround(CELLS, vert))
			for i in reasonablyaround(CELLS,xn)
				if i!=vert && distance(xn, surf.vertices[Int(i)])<2*site_radius
					good=false
				end
			end
		end
		if good
			DE, surf=ev_envm(surf, vert, xn, k, part, mu) #here is done also check for maximum dihedral angle
			p=minimum([1., exp(-b*DE)])
			if rand() <p
				#accept
				#surf.vertices[vert]=xn
				#println("moved vertex ", vert)
				#surf = update_geoquantvm(surf, vert, xn)
				CELLS = update_cellsvm(CELLS, vert, xn)
				acc+=1
			else
				#reject
				surf = update_geoquantvm(surf, vert, x0)
			end
		end
	end
	if verbose
		println("acceptability for vertex move: " , acc/Nsvert)
	end
	a[1]=acc/Nsvert
	#mc step for links
	acc=0
	for steplink in 1:Nslink
		i=Int16(rand(1:length(surf.edges)))
		t1,t2=surf.edge_faces[i][1],surf.edge_faces[i][2] #triangles (faces) sharing edge i
		e1,e2=surf.edges[i][1],surf.edges[i][2]
		ne1=filter(x -> (x!= e1 && x!= e2), surf.faces[t1])[1]
		ne2=filter(x -> (x!= e1 && x!= e2), surf.faces[t2])[1]
		good=true
		#first check if flipping this edge still satisfies maximum edge length constraints
		if distance(surf.vertices[ne1],surf.vertices[ne2]) > max_edge_len
			good=false
		end
		#check for pyramidal structures and edge conservation will be done inside update_geoquantlm() funciton
		if good	
			DE, surf, check=ev_enlm(surf, i, k, part, mu, [e1,e2,ne1,ne2]) #here also check for pyramidal structures and edge conservation
			if check
				p=minimum([1., exp(-b*DE)])
				if rand()<p
					#accept
					acc+=1
				else
					#reject and restore previous configuration
					surf, = update_geoquantlm(surf, i, [e1,e2,ne1,ne2])
				end
			else
				#reject without doing nothing
			end
		end
	end
	if verbose
		println("acceptability for link move: ", acc/Nslink)
	end
	a[2]=acc/Nslink
	acc=0
	#MC sweep for particles
	for steppart in 1:Nspart
		v=rand(1:surf.Nv)
		ov=deepcopy(part[v])
		goto=rand(surf.neig[v])
		ogoto=deepcopy(part[goto])
		if ov!=ogoto #if they're equal, DeltaH=0 => accept always
			Ein=ev_energyonset(surf,k,part,mu,[v,goto])
			part[goto]=ov
			part[v]=ogoto
			Efin=ev_energyonset(surf,k,part,mu,[v,goto])
			DE=Efin-Ein
			p=minimum([1., exp(-b*DE)])
			if rand() < p
				#accept
				acc+=1
			else
				#reject, undo particle switch
				part[goto]=ogoto
				part[v]=ov
			end
		end
	end
	if verbose
		println("acceptability for particle move (only nontrivial exchange counted): ", acc/Nspart)
	end
	a[3]=acc/Nspart
	return surf, part, a
end

function NsweepsPM(surf::SurfGeoQuant, CELLS, sigma, max_edge_len, site_radius, k, part, mu, Ki, Kd, g, N, filename="default", M=0,savecfg=false)
	for i in 1:N
		println("step ", i)
		@time surf, part=MC_sweepPM(surf,CELLS,sigma,max_edge_len,site_radius,k,part,mu, Ki, Kd, g)
		if savecfg==true && i%M==0
			save_config(surf, part, filename*"withparticles_step"*string(i))
		end
	end
	return surf, part
end
