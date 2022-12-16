using Random
Random.seed!(1234)
include("sphere.jl") #requires TriSurface and CellsGrid definition
I3=[1. 0. 0.; 0. 1. 0.; 0. 0. 1.] #identity 3x3


function ev_curv(su::SurfGeoQuant,n)
	Av=0
	Nv=[0.,0.,0.]
	for nf in su.neig_faces[n]
		af=su.face_area[nf]
		Av+=af
		Nv=sum_vec(Nv,(su.face_normals[nf]*af))
	end
	Nv/=norm(Nv)
	Av/=3
	Sv=zeros(3,3)
	Pv=I3-tensor_prod(Nv,Nv)
	for e in su.neig_edges[n]
		#su.edges[e][1] == n ? j=su.edges[e][2] : j=su.edges[e][1] # j = neighbour of i, with which edge e is in common
		Ne=sum_vec(su.face_normals[su.edge_faces[e][2]],su.face_normals[su.edge_faces[e][1]])
		Ne/=norm(Ne)
		We=scalar_prod(Nv, Ne)
		Sv=weighted_sum(Sv,We,su.Se[e])
	end
	Sv=Pv*Sv*Pv #notice that transpose(Pv)=Pv
	Sv/=Av
	cm=(Sv[1,1]+Sv[2,2]+Sv[3,3])/2
	lc=0 #local curvature
	if part[n]==+1 #if a particle is present on vertex, then induce local curvature
		lc=mu
	end
	return cm
end

function ev_energyonset(su::SurfGeoQuant, k, part, mu, set)
	E=0
	for n in set
		Av=0
		Nv=[0.,0.,0.]
		for nf in su.neig_faces[n]
			af=su.face_area[nf]
			Av+=af
			Nv=sum_vec(Nv,(su.face_normals[nf]*af))
		end
		Nv/=norm(Nv)
		Av/=3
		#Sv=zeros(3,3)
		Sv=Float64[0. 0. 0.;0. 0. 0.;0. 0. 0.]
		Pv=I3-tensor_prod(Nv,Nv)
		for e in su.neig_edges[n]
			#su.edges[e][1] == n ? j=su.edges[e][2] : j=su.edges[e][1] # j = neighbour of i, with which edge e is in common
			Ne=sum_vec(su.face_normals[su.edge_faces[e][2]],su.face_normals[su.edge_faces[e][1]])
			Ne/=norm(Ne)
			We=scalar_prod(Nv, Ne)
			Sv=weighted_sum(Sv,We,su.Se[e])
		end
		Sv=Pv*Sv*Pv #notice that transpose(Pv)=Pv and Pv is made of real number => adjoint(Pv)=Pv
		#Householder transformation used to evaluate Sv eigenvalues faster (test against LinearAlgebra.eigvals, is much faster
		#c1,c2=eigenHH(Sv, Nv)
		#println(c1+c2)
		#println(Sv[1,1]+Sv[2,2]+Sv[3,3])
		#cm=(c1+c2)/2
		cm=(Sv[1,1]+Sv[2,2]+Sv[3,3])/(2*Av)
		lc=0 #local curvature
		if part[n]==+1 #if a particle is present on vertex, then induce local curvature
			lc=mu
		end
		E+=((cm-lc)^2)*Av
	end
	return E*k*0.5
end

#function to evaluate total energy, OLD NEEDS TO BE FIXED
function ev_energy(su::SurfGeoQuant, k, part, mu)
	En=ev_energyonset(su,k,part,mu, collect(1:length(su.vertices)))
	return En
end


#evaluate energy difference in case of moving vertex i to new position xn (assumin this is an acceptable position)
function ev_envm(su::SurfGeoQuant, i, xn, k, part, mu) 
	#when vertex i is displaced to evaluate energy difference we need to evaluate change in curvatures for i and all his neighbours
	Ein=ev_energyonset(su, k, part, mu, vcat(su.neig[i], [i]))
	su2 = update_geoquantvm(su, i, xn)
	Efin=ev_energyonset(su2,k,part,mu, vcat(su2.neig[i], [i]))
	return Efin-Ein, su2
end

function ev_enlm(su::SurfGeoQuant, e, k, part, mu)
	t1,t2=su.edge_faces[e][1],su.edge_faces[e][2]
	inv_ver=[su.faces[t1]; su.faces[t2]] #vertices involved whose energy will change, to concatenate two vectors a and b used [a; b]
	inv_ver=union(inv_ver)
	Ein=ev_energyonset(su,k,part,mu,inv_ver)
	su2, g = update_geoquantlm(su, e)
	if g==false
		return 0, su2, g
	end
	Efin=ev_energyonset(su2,k,part,mu,inv_ver)
	return Efin-Ein, su2, g
end

#function MC sweep accounting for particles move
function MC_sweepPM(surf::SurfGeoQuant, CELLS, sigma, max_edge_len, site_radius, k, part, mu,Nsvert, Nslink, Nspart; verbose=false)
	b=1.
	#mc step for every vertex
	acc=0
	for stepvertex in 1:Nsvert
		vert=rand(collect(1:surf.Nv))
		x0=deepcopy(surf.vertices[vert])
		xn=sum_vec(x0,((rand(3).-0.5)*2*sigma))
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
	#mc step for links
	acc=0
	for steplink in 1:Nslink
		i=rand(1:length(surf.edges))
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
			DE, surf, check=ev_enlm(surf, i, k, part, mu) #here also check for pyramidal structures and edge conservation
			if check
				p=minimum([1., exp(-b*DE)])
				if rand()<p
					#accept
					acc+=1
				else
					#reject and restore previous configuration
					surf, = update_geoquantlm(surf, i)
				end
			else
				#reject without doing nothing
			end
		end
	end
	if verbose
		println("acceptability for link move: ", acc/Nslink)
	end
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
				println("accepted particle move with DE=", DE)
				acc+=1
			else
				#reject, undo particle switch
				part[goto]=ogoto
				part[v]=ov
			end
		end
	end
	if verbose
		println("acceptability for particle move: ", acc/Nspart)
	end
	return surf, part
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

