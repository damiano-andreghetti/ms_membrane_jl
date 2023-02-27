using Random
include("sphere.jl") #requires TriSurface and CellsGrid definition
function ev_energyonset(su::SurfGeoQuant,Ene::Vector{Float64}, k::Float64, sigma::Float64, part::Vector{Bool}, A::Float64, B::Float64,set::Vector{Int})
	for n in set
		cm, Av =ev_curv_vertex(su, n)
		Ene[n]=(k*cm*cm + sigma)*Av
		if part[n]
			Ene[n]+=((A*cm+B*cm*cm)*Av)
		end
	end
	return Ene
end

#function to evaluate total energy, OLD NEEDS TO BE FIXED
function ev_energy(su::SurfGeoQuant, k::Float64, sigma::Float64,part::Vector{Bool}, A::Float64,B::Float64)
	En=ev_energyonset(su,zeros(su.Nv),k*0.5,sigma,part,A,B, [Int(i) for i in 1:length(su.vertices)])
	return En
end


#evaluate energy difference in case of moving vertex i to new position xn (assumin this is an acceptable position)
function ev_envm(su::SurfGeoQuant, Ene::Vector{Float64}, i::Int, xn::MVector{3,Float64}, k::Float64, sigma::Float64,part::Vector{Bool}, A::Float64,B::Float64) 
	#when vertex i is displaced to evaluate energy difference we need to evaluate change in curvatures for i and all his neighbours
	su = update_geoquantvm(su, i, xn)
	Ene=ev_energyonset(su,Ene,k,sigma,part,A, B,vcat([i],su.neig[i]))
	return su,Ene
end

function ev_enlm(su::SurfGeoQuant, Ene::Vector{Float64}, e::Int, k::Float64, sigma::Float64,part::Vector{Bool}, A::Float64, B::Float64,set::Vector{Int})
	su, g = update_geoquantlm(su, e, set)
	if g==false
		return su, g, Ene
	end
	Ene=ev_energyonset(su,Ene,k,sigma,part,A,B,set)
	return su, g, Ene
end

#function MC sweep accounting for particles move
function MC_sweepPM(surf::SurfGeoQuant, CELLS::CellsGrid, Ene::Vector{Float64},l0::Float64, max_edge_len::Float64, site_radius::Float64, k::Float64, sigma::Float64, part::Vector{Bool}, A::Float64,B::Float64,Nsvert::Int64, Nslink::Int64; verbose::Bool=false)
	b=1.
	k=k*0.5 #from now on we'll just need k/2 so there's no matter of repeating this operation many times
	#mc step for every vertex
	acc=0
	cont=0
	for stepvertex in 1:Nsvert
		vert=Int(rand(1:surf.Nv))
		x0=deepcopy(surf.vertices[vert])
		#xn=sum_vec(x0,((rand(3).-0.5)*2*sigma))
		delta=MVector{3,Float64}((rand(3).-0.5)*2*l0)
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
			cont+=1
			Ein=sum(Ene[surf.neig[vert]])
			Ein+=Ene[vert]
			surf, Ene=ev_envm(surf, Ene, vert, xn, k, sigma,part, A, B) #here is done also check for maximum dihedral angle
			DE=sum(Ene[surf.neig[vert]])+Ene[vert]-Ein
			p=minimum([1., exp(-b*DE)])
			if rand() <p
				#accept
				CELLS = update_cellsvm(CELLS, vert, xn)
				acc+=1
			else
				#reject
				surf, Ene=ev_envm(surf,Ene,vert,x0,k,sigma,part,A,B)
			end
		end
	end
	if verbose
		println("acceptance for vertex moves (considering only feasible moves): " , acc/cont)
	end
	#mc step for links
	acc=0
	cont=0
	for steplink in 1:Nslink
		i=Int(rand(1:length(surf.edges)))
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
			Ein=Ene[e1]+Ene[e2]+Ene[ne1]+Ene[ne2]
			surf, check, Ene=ev_enlm(surf, Ene, i, k, sigma,part, A, B,[e1,e2,ne1,ne2]) #here also check for pyramidal structures and edge conservation
			if check
				cont+=1
				DE=Ene[e1]+Ene[e2]+Ene[ne1]+Ene[ne2]-Ein
				p=minimum([1., exp(-b*DE)])
				if rand()<p
					#accept
					acc+=1
				else
					#reject and restore previous configuration
					surf,dumbcheck,Ene = ev_enlm(surf, Ene, i, k,sigma, part, A, B,[e1,e2,ne1,ne2])
				end
			else
				#reject without doing nothing
			end
		end
	end
	if verbose
		println("acceptance for link moves (considering only feasible moves): ", acc/cont)
	end
	return surf, part, Ene
end
