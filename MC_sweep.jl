using Random
Random.seed!(1234)
include("sphere.jl") #requires TriSurface and CellsGrid definition
I3=[1. 0. 0.; 0. 1. 0.; 0. 0. 1.] #identity 3x3

function eigenCv(C)
	p=C[2,2]
	q=C[2,3]
	r=C[3,2]
	s=C[3,3]
	c1=0
	c2=0
	if r!=0 && q!=0
		d=(s+p)^2-4*(p*s-r*q)
		if d>0
			c1=0.5*((s+p)+sqrt(d))
			c2=0.5*((s+p)-sqrt(d))
		else
			c1=c2=(p+s)*0.5
		end
	else
		if p>s
			c1=p
			c2=s
		else
			c1=s
			c2=p
		end
	end
	return c1, c2
end

function eigenHH(Sv, Nv)
	x=[1.,0.,0.]
	wp=sum_vec(x,Nv)
	wm=sub_vec(x,Nv)
	nwp=norm(wp)
	nwm=norm(wm)
	if nwp>=nwm
		wp/=nwp
		H=weighted_sum(I3,-2.,tensor_prod(wp,wp))
	else
		wm/=nwm
		H=weighted_sum(I3,-2.,tensor_prod(wm,wm))
	end
	Cv=transpose(H)*Sv*H
	c1,c2 = eigenCv(Cv)
	return c1, c2
end

function ev_energyonset(su::SurfGeoQuant, k, part, mu, set)
	E=0
	if 1 in set
		open("myfile.txt", "a+") do io
			println(io, "evaluating energy on a set of vertices with a particle, set=", set)
		end
	end
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
		Sv=zeros(3,3)
		Pv=I3-tensor_prod(Nv,Nv)
		for e in su.neig_edges[n]
			#su.edges[e][1] == n ? j=su.edges[e][2] : j=su.edges[e][1] # j = neighbour of i, with which edge e is in common
			Ne=sum_vec(su.face_normals[su.edge_faces[e][2]],su.face_normals[su.edge_faces[e][1]])
			Ne/=norm(Ne)
			We=scalar_prod(Nv, Ne)
			Sv=weighted_sum(Sv,We,su.Se[e])
		end
		Sv=transpose(Pv)*Sv*Pv
		Sv/=Av
		#Householder transformation used to evaluate Sv eigenvalues faster (test against LinearAlgebra.eigvals, is much faster
		#c1,c2=eigenHH(Sv, Nv)
		#println(c1+c2)
		#println(Sv[1,1]+Sv[2,2]+Sv[3,3])
		#cm=(c1+c2)/2
		cm=(Sv[1,1]+Sv[2,2]+Sv[3,3])/2
		lc=0 #local curvature
		if part[n]==+1 #if a particle is present on vertex, then induce local curvature
			lc=mu
			open("myfile.txt", "a+") do io
				println(io, "particle present, Hc=",cm, " H0=",lc," Hc-H0=",cm-lc," (Hc-H0)^2=", (cm-lc)^2, "energy for this vertex is=",(k)*0.5*((cm-lc)^2)*Av)
			end
		end
		E+=(k)*0.5*((cm-lc)^2)*Av
	end
	return E
end

#function to evaluate total energy, OLD NEEDS TO BE FIXED
function ev_energy(su::SurfGeoQuant, k, part, mu)
	En=ev_energyonset(su,k,part,mu, collect(1:length(su.vertices)))
	return En
end


#evaluate energy difference in case of moving vertex i to new position xn (assumin this is an acceptable position)
function ev_envm(su::SurfGeoQuant, i, xn, k, part, mu, maxang) 
	#when vertex i is displaced to evaluate energy difference we need to evaluate change in curvatures for i and all his neighbours
	if part[i]==1
		open("myfile.txt", "a+") do io
			println(io, "evaluating energy change for a move of a vertex with a particle")
		end
	end
	Ein=ev_energyonset(su, k, part, mu, vcat(su.neig[i], [i]))
	su2, g = update_geoquantvm(su, i, xn, maxang)
	if g==false
		open("myfile.txt", "a+") do io
			println(io, "dihedral angle constraint violated, vertex move interrupted")
		end
		return 0, su2, g
	end
	Efin=ev_energyonset(su2,k,part,mu, vcat(su2.neig[i], [i]))
	return Efin-Ein, su2, g
end

function ev_enlm(su::SurfGeoQuant, e, k, part, mu, maxang)
	t1,t2=su.edge_faces[e][1],su.edge_faces[e][2]
	inv_ver=[su.faces[t1]; su.faces[t2]] #vertices involved whose energy will change, to concatenate two vectors a and b used [a; b]
	union(inv_ver)
	Ein=ev_energyonset(su,k,part,mu,inv_ver)
	su2, g = update_geoquantlm(su, e, maxang)
	if g==false
		return 0, su2, g
	end
	Efin=ev_energyonset(su2,k,part,mu,inv_ver)
	return Efin-Ein, su2, g
end

#function MC sweep accounting for particles move
function MC_sweepPM(surf::SurfGeoQuant, CELLS, sigma, max_edge_len, site_radius, k, part, mu, maxang,Nsvert, Nslink, Nspart; verbose=false)
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
			DE, surf, angcheck=ev_envm(surf, vert, xn, k, part, mu, maxang) #here is done also check for maximum dihedral angle
			p=minimum([1., exp(-b*DE)])
			if rand() <p && angcheck 
				#accept
				#surf.vertices[vert]=xn
				#println("moved vertex ", vert)
				#surf = update_geoquantvm(surf, vert, xn)
				CELLS = update_cellsvm(CELLS, vert, xn)
				acc+=1
			else
				#reject
				surf,  = update_geoquantvm(surf, vert, x0, maxang)
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
		#first check if flipping this edge still satisfies non overlapping and maximum edge length constraints
		if 2*site_radius > distance(surf.vertices[ne1],surf.vertices[ne2]) || distance(surf.vertices[ne1],surf.vertices[ne2]) > max_edge_len
			good=false
			#println(distance(surf.vertices[ne1],surf.vertices[ne2]))
		end
		if length(surf.neig_edges[surf.edges[i][1]])-1<3 || length(surf.neig_edges[surf.edges[i][2]])-1<3
			#println("doesnt satisfy edge conservation")
			good=false
		end
		if ne1 in surf.neig[ne2]
			#if edge already exist dont accept move (pyramidal structures)
			good=false
		end
		if good
			DE, surf, angcheck=ev_enlm(surf, i, k, part, mu, maxang) #here also check for maximum dihedral angle constraint
			p=minimum([1., exp(-b*DE)])
			if angcheck && rand() < p
				#accept
				#surf = update_geoquantlm(surf, i)
				acc+=1
			else
				#reject
				surf, = update_geoquantlm(surf, i, maxang)
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
		if ov!=ogoto #if they're equal, DeltaH=0 => tanh(0)=0 => reject always
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

