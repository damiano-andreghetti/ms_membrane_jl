using Random
include("sphere.jl") #requires TriSurface and CellsGrid definition
I3=[1 0 0; 0 1 0; 0 0 1] #identity 3x3

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
	x=[1,0,0]
	wp=x.+Nv
	wm=x.-Nv
	if norm(wp)>=norm(wm)
		wp/=norm(wp)
		H=I3-2*tensor_prod(wp,wp)
	else
		wm/=norm(wm)
		H=I3-2*tensor_prod(wm,wm)
	end
	Cv=adjoint(H)*Sv*H
	c1,c2 = eigenCv(Cv)
	return c1, c2
end

function ev_energyonset(su::SurfGeoQuant, k, part, mu, set)
	E=0
	for n in set
		Av=0
		Nv=[0.,0.,0.]
		for nf in su.neig_faces[n]
			af=su.face_area[nf]
			Av+=af
			Nv=Nv.+(su.face_normals[nf]*af)
		end
		Nv/=norm(Nv)
		Av/=3
		Sv=zeros(3,3)
		Pv=I3-tensor_prod(Nv,Nv)
		for e in su.neig_edges[n]
			su.edges[e][1] == n ? j=su.edges[e][2] : j=su.edges[e][1] # j = neighbour of i, with which edge e is in common
			ed1=[] #TO FIND FACES IN COUNTERCLOCKWISE ORDER, S.T. DIHEDRAL ANGLE SIGN IS CORRECT (give cross product involved in a consistent order)
            ed2=[]
            for z in 1:3 #uses the fact that triangle indices are spatially in a counterclockwise order
                if su.faces[su.edge_faces[e][1]][z]==n #ed will be the two "border" vertices of the hexagon in counterclockwise order
                    push!(ed1, su.faces[su.edge_faces[e][1]][((z+1)%3)+1])
					push!(ed1, su.faces[su.edge_faces[e][1]][((z+2)%3)+1])
				end
                if su.faces[su.edge_faces[e][2]][z]==n
                    push!(ed2, su.faces[su.edge_faces[e][2]][((z+1)%3)+1])                    
                    push!(ed2, su.faces[su.edge_faces[e][2]][((z+2)%3)+1])
				end
			end
			if length(ed2)<2
				println(ed1, ed2)
				println(su.faces[su.edge_faces[e][2]])
				println(e)
				println(su.edges[e])
			end
            if ed1[2]==ed2[1] #if they are already in counter clockwise order
                Nf_1=su.face_normals[su.edge_faces[e][1]]
                Nf_2=su.face_normals[su.edge_faces[e][2]]
            else #if not exchange them
				Nf_1=su.face_normals[su.edge_faces[e][2]]
                Nf_2=su.face_normals[su.edge_faces[e][1]]
			end
			ae=0
			if part[su.edge_faces[e][1]]==1
				ae+=1
			end
			if part[su.edge_faces[e][2]]==1
				ae+=1
			end
			Ne=Nf_1.+Nf_2
			Ne/=norm(Ne)
			re=su.vertices[n].-su.vertices[j]
			be=cross_prod(re/norm(re),Ne)
			phi=sign(sum(cross_prod(Nf_1, Nf_2).*re))*acos(sum(Nf_1.*Nf_2))+pi 
			he=2*norm(re)*cos(phi/2)
			Se=(he-ae*mu)*tensor_prod(be,be) #MODIFY here to account for particles
			We=sum(Nv.*Ne) #dot produtc between Nv and Ne
			Sv+=(We*adjoint(Pv)*Se*Pv)
		end
		Sv/=Av
		#Householder transformation used to evaluate Sv eigenvalues faster (test against LinearAlgebra.eigvals, is much faster
		c1,c2=eigenHH(Sv, Nv)
		E+=k*0.5*((c1+c2)^2)*Av
	end
	return E
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
	return Efin-Ein
end

function ev_enlm(su::SurfGeoQuant, e, k, part, mu)
	t1,t2=su.edge_faces[e][1],su.edge_faces[e][2]
	inv_ver=[su.faces[t1]; su.faces[t2]] #vertices involved whose energy will change, to concatenate two vectors a and b used [a; b]
	union(inv_ver)
	Ein=ev_energyonset(su,k,part,mu,inv_ver)
	su2 = update_geoquantlm(su, e)
	Efin=ev_energyonset(su2,k,part,mu,inv_ver)
	return Efin-Ein
end

function MC_sweep(surf::SurfGeoQuant, CELLS, sigma, max_edge_len, site_radius, k, part, mu, b)
	#mc step for every vertex
	vert_order=shuffle(collect(1:surf.Nv))
	acc=0
	for vert in vert_order
		x0=copy(surf.vertices[vert])
		xn=x0.+((rand(3).-0.5)*2*sigma)
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
		if good
			for i in reasonablyaround(CELLS, vert)
				#println(distance(xn, surf.vertices[Int(i)]))
				if good==false
					break;
				end
				if distance(xn, surf.vertices[Int(i)])<2*site_radius
					good=false
				end
			end
		end
		if good
			DE=ev_envm(surf, vert, xn, k, part, mu)
			p=minimum([1., exp(-b*DE)])
			if rand() < p
				#accept
				#surf.vertices[vert]=xn
				#println("moved vertex ", vert)
				surf = update_geoquantvm(surf, vert, xn)
				CELLS = update_cellsvm(CELLS, vert, xn)
				acc+=1
			else
				#reject
			end
		end
	end
	println("acceptability for vertex move: " , acc/length(vert_order))
	#mc step for links
	acc=0
	edge_order=shuffle(collect(1:length(surf.edges)))[1:surf.Nv]
	for i in edge_order
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
		if good
			DE=ev_enlm(surf, i, k, part, mu)
			p=minimum([1., exp(-b*DE)])
			if rand() < p
				#accept
				surf = update_geoquantlm(surf, i)
				acc+=1
			else
				#reject
			end
		end
	end
	println("acceptability for link move: ", acc/length(edge_order))
	return surf
end

function Nsweeps(surf::SurfGeoQuant, CELLS, sigma, max_edge_len, site_radius, k, part, mu, N, beta)
	for i in 1:N
		@time surf=MC_sweep(surf,CELLS,sigma,max_edge_len,site_radius,k,part,mu, beta)
		println(ev_energy(surf,k,part,mu))
	end
	return surf
end
