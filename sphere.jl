using Meshes #requires Meshes package for perp4plot function
#various function used for vector operations
function norm(a::Vector{})
	return sqrt(a[1]^2+a[2]^2+a[3]^2)
end
function normalize_vec(v)
	return v/norm(v)
end
function normalize_vec(v, L)
	return v*L/norm(v)
end
function mpf(v1,v2)
	nv=((v1.+v2)*0.5)
	return nv/norm(nv)
end
function distance(a,b)
	return norm(a.-b)
end	
function cross_prod(a,b)
	return [a[2]*b[3]-a[3]b[2], a[3]b[1]-a[1]b[3], a[1]b[2]-a[2]b[1]]
end
function area(AB, AC)
	return norm(cross_prod(AB, AC))*0.5
end
#tensor product between two vectors 
#actually corresponds to outer product: results in a matrix A s.t. Aij=a_i*b_j
function tensor_prod(a,b)
	A=zeros(3,3)
	for i in 1:3
		for j in 1:3
			A[i,j]=a[i]*b[j]
		end
	end
	return A
end
#function to check in faces are counterclockwise
#find center of each triangular face[A,B,C], call it O (O represent vector from origin to O)
#evaluate cross product v=AB X AC, then dot(v,O)
#if sign of this dot productt is positive then vertices are counter clockwise
function chechkcounterclock(vertices)
	A=vertices[1]
	B=vertices[2]
	C=vertices[3]
	O=A.+B
	O=O.+C
	O/=3
	AB=B.-A
	AC=C.-A
	vect=cross_prod(AB,AC)
	return Bool(sum(vect.*O)>0)
end
#function used to evaluate size of variables (for newly defined types)
#REMOVED: use Base.summarysize(x) instead


#
#
#TriSurface definition and related functions
#
#
#struct for triangulated surface, just surface definition, no related geometrical quantities
struct TriSurface
	Nv::Int
	vertices::Vector{Vector{Number}}
	faces::Vector{Vector{Number}}
	A::BitArray{}
	neig::Vector{Vector{Number}}
end
#subdivide a given trisurface: each triangle is divided into 4 triangle
function subdivide(surf::TriSurface)
	newvertices = copy(surf.vertices)
	edges=[]
	for el in surf.faces
		sg1=[el[1],el[2]]
		sg2=[el[2],el[3]]
		sg3=[el[1],el[3]]
		sg4=[el[2],el[1]]
		sg5=[el[3],el[2]]
		sg6=[el[3],el[1]]
		push!(edges, sg1)
		push!(edges, sg2)
		push!(edges, sg3)
		push!(edges, sg4)
		push!(edges, sg5)
		push!(edges, sg6)
	end
	edges = filter(x->x[1]<x[2],edges)
	edges=union(edges)
	epoint(v1,v2) = findfirst(x->x==(v2>v1 ? [v1,v2] : [v2,v1]), edges) + length(surf.vertices)

	newfaces = []
	resize!(newvertices,length(surf.vertices) + length(edges))
	#new faces added keep count-clockwise vertex order (only if they originally had)
	for f in surf.faces
		v1=f[1]
		v2=f[2]
		v3=f[3]
		ev3 = epoint(v1,v2)
		ev1 = epoint(v2,v3)
		ev2 = epoint(v3,v1)

		newvertices[ev3] = mpf(surf.vertices[v1],surf.vertices[v2])
		newvertices[ev1] = mpf(surf.vertices[v2],surf.vertices[v3])
		newvertices[ev2] = mpf(surf.vertices[v3],surf.vertices[v1])

		push!(newfaces, [v1,ev3,ev2])
		push!(newfaces, [v2,ev1,ev3])
		push!(newfaces, [v3,ev2,ev1])
		push!(newfaces, [ev1,ev2,ev3])
	end
	edges=[]
	for el in newfaces
		sg1=[el[1],el[2]]
		sg2=[el[2],el[3]]
		sg3=[el[1],el[3]]
		sg4=[el[2],el[1]]
		sg5=[el[3],el[2]]
		sg6=[el[3],el[1]]
		push!(edges, sg1)
		push!(edges, sg2)
		push!(edges, sg3)
		push!(edges, sg4)
		push!(edges, sg5)
		push!(edges, sg6)
	end
	edges = filter(x->x[1]<x[2],edges)
	edges=union(edges)
	Nv=length(newvertices)
	A=zeros(Bool, Nv, Nv)
	for e in edges
		A[e[1],e[2]]=1
		A[e[2],e[1]]=1
	end
	neig=[]
	for i in 1:Nv
		push!(neig, findall(x->x==true,A[i,:]))
	end
	return TriSurface(Nv,newvertices, newfaces, A, neig)
end
#iterated call for previous function
function subNtimes(Nsub,surf::TriSurface)
	for i in 1:Nsub
		surf= subdivide(surf)
	end
	return surf
end
#creates an icosahedron: 12 vertices, 20 faces,  is used to create initial pseudo-spherical mesh using then subdivide function 
function create_ico()
	t = ( 1 + sqrt( 5 ) ) / 2; #sezione aurea
	ico_vertices = [[ -1,  t,  0 ], [  1, t, 0 ], [ -1, -t,  0 ], [  1, -t,  0 ],
					[  0, -1,  t ], [  0, 1, t ], [  0, -1, -t ], [  0,  1, -t ],
					[  t,  0, -1 ], [  t, 0, 1 ], [ -t,  0, -1 ], [ -t,  0,  1 ]]

	ico_faces = [[1, 12, 6], [1, 6, 2], [1, 2, 8], [1, 8, 11], [1, 11, 12], [2, 6, 10], [6, 12, 5],
				[12, 11, 3], [11, 8, 7], [8, 2, 9], [4, 10, 5], [4, 5, 3], [4, 3, 7], [4, 7, 9],
				[4, 9, 10], [5, 10, 6], [3, 5, 12], [7, 3, 11], [9, 7, 8], [10, 9, 2]]
	
	ico_vertices= normalize_vec.(ico_vertices)
	
	edges=[]
	for el in ico_faces
		sg1=[el[1],el[2]]
		sg2=[el[2],el[3]]
		sg3=[el[1],el[3]]
		sg4=[el[2],el[1]]
		sg5=[el[3],el[2]]
		sg6=[el[3],el[1]]
		push!(edges, sg1)
		push!(edges, sg2)
		push!(edges, sg3)
		push!(edges, sg4)
		push!(edges, sg5)
		push!(edges, sg6)
	end
	edges = filter(x->x[1]<x[2],edges)
	edges=union(edges)
	Nv = length(ico_vertices)
	A=zeros(Bool, Nv,Nv)
	for e in edges
		A[e[1],e[2]]=1
		A[e[2],e[1]]=1
	end
	neig=[]
	for i in 1:Nv
		push!(neig, findall(x->x==true,A[i,:]))
	end
	return TriSurface(Nv, ico_vertices, ico_faces, A, neig)
end
#used to adjust initial average edge length, this function sets it equal to 1.2, being able to have all initial edges between 1 and sqrt(2)
function normEdges(surf::TriSurface)
	av=0
	edges=[]
	for el in surf.faces
		sg1=[el[1],el[2]]
		sg2=[el[2],el[3]]
		sg3=[el[1],el[3]]
		sg4=[el[2],el[1]]
		sg5=[el[3],el[2]]
		sg6=[el[3],el[1]]
		push!(edges, sg1)
		push!(edges, sg2)
		push!(edges, sg3)
		push!(edges, sg4)
		push!(edges, sg5)
		push!(edges, sg6)
	end
	edges = filter(x->x[1]<x[2],edges)
	edges=union(edges)
	lmin=1
	lmax=1
	for e in edges
		av+=distance(surf.vertices[e[1]],surf.vertices[e[2]])
		if distance(surf.vertices[e[1]],surf.vertices[e[2]]) >lmax
			lmax=distance(surf.vertices[e[1]],surf.vertices[e[2]])
		end
		if distance(surf.vertices[e[1]],surf.vertices[e[2]]) <lmin
			lmin=distance(surf.vertices[e[1]],surf.vertices[e[2]])
		end
	end
	av/=length(edges)
	println("average distance between vertices is: ", av)
	println("minimum is ", lmin, " maiximum is ", lmax)
	c=1.2/av
	surf.vertices.*=c
	av=0
	lmin=1
	lmax=1
	for e in edges
		av+=distance(surf.vertices[e[1]],surf.vertices[e[2]])
		if distance(surf.vertices[e[1]],surf.vertices[e[2]]) >lmax
			lmax=distance(surf.vertices[e[1]],surf.vertices[e[2]])
		end
		if distance(surf.vertices[e[1]],surf.vertices[e[2]]) <lmin
			lmin=distance(surf.vertices[e[1]],surf.vertices[e[2]])
		end
	end
	av/=length(edges)
	println("now is: ", av)
	println("minimum is ", lmin, " maiximum is ", lmax)
	return TriSurface(surf.Nv, surf.vertices, surf.faces, surf.A, surf.neig)
end


#
#
#CellsGrid definition and related functions
#
#
#struct for dividing simulation space in cells (to avoid O(N^2) MC sweep)
struct CellsGrid
	L::Number
	Nx::Number
	Ny::Number
	Nz::Number
	cells::Vector{Vector{Vector{Vector{Number}}}}
	vertcell::Vector{Vector{Number}}
end
#function to generate cells, given a cube of side L centered in origin, it is divided into Nx by Ny by Nz cells
function genCells(surf::TriSurface, Nx,Ny,Nz, L)
	ver=reduce(vcat,transpose.(surf.vertices))
	Lx=L/Nx
	Ly=L/Ny
	Lz=L/Nz
	vertcell=[[0,0,0] for i in 1:(length(surf.vertices))] #vector of pointers for each particle to preceeding particle in same cell. -1 if no other particle
	cells=[[[[] for k in 1:Nz] for j in 1:Ny] for i in 1:Nx] #vector of headers containing last particle inserted in each cell
	ver=surf.vertices
	for i in 1:length(ver)
		cx, cy, cz= Int(floor((ver[i][1]+L*0.5)/Lx) + 1), Int(floor((ver[i][2]+L*0.5)/Ly) + 1), Int(floor((ver[i][3]+L*0.5)/Lz)+ 1)#floor division to find cartesian coordinates of cell
		#println(cx,cy,cz)
		push!(cells[cx][cy][cz], i)
		vertcell[i]=[cx,cy,cz]
	end
	return CellsGrid(L, Nx,Ny,Nz,cells,vertcell)
end
#update cells after vertex move
function update_cellsvm(ce::CellsGrid, i, xn)
	filter!(e->e==i,ce.cells[ce.vertcell[i][1]][ce.vertcell[i][2]][ce.vertcell[i][3]])
	L=ce.L
	Nx=ce.Nx
	Ny=ce.Ny
	Nz=ce.Nz
	Lx=L/Nx
	Ly=L/Ny
	Lz=L/Nz
	cx, cy, cz= Int(floor((xn[1]+L*0.5)/Lx) + 1), Int(floor((xn[2]+L*0.5)/Ly) + 1), Int(floor((xn[3]+L*0.5)/Lz)+ 1)
	push!(ce.cells[cx][cy][cz], i)
	ce.vertcell[i]=[cx,cy,cz]
	return ce
end
#returns vertices reasonably around i to check for overlapping
function reasonablyaround(ce::CellsGrid, i)
	ix=ce.vertcell[i][1]
	iy=ce.vertcell[i][2]
	iz=ce.vertcell[i][3]
	rsn=ce.cells[ix][iy][iz]
	#now include also neighbour cells
	#use pbc for this
	for vx in [-1,0,1]
		for vy in [-1,0,1]
			for vz in [-1, 0, 1]
				nx=ix+vx
				ny=iy+vy
				nz=iz+vz
				if (nx>=1 && nx <=ce.Nx) && (ny>=1 && ny <=ce.Ny) && (nz>=1 && nz <=ce.Nz)
					vcat(rsn, ce.cells[nx][ny][nz])
				end
			end
		end
	end
	rsn=collect(1:length(ce.vertcell))
	rsn=union(rsn)
	filter!(e->e≠i,rsn) #remove vertex i from this list
	return rsn
end

#
#
#SurfGeoQuant definition and related functions
#
#
#struct for triangulated surface, but equipped with some geometrical quantities that will result useful for elastic energy calculation
struct SurfGeoQuant
	Nv::Int
	vertices::Vector{Vector{Number}}
	faces::Vector{Vector{Number}}
	edges::Vector{Vector{Number}}
	neig::Vector{Vector{Number}}
	neig_faces::Vector{Vector{Number}} #faces around each vertex
	neig_edges::Vector{Vector{Number}} #edges around each vertex
	face_area::Vector{Number} #area of faces
	edge_faces::Vector{Vector{Number}} #faces around each edge
	face_normals::Vector{Vector{Number}} #face normals
end

#define copy for surfgeoquant
Base.deepcopy(su::SurfGeoQuant) = SurfGeoQuant(deepcopy(su.Nv), deepcopy(su.vertices), deepcopy(su.faces), deepcopy(su.edges), deepcopy(su.neig), deepcopy(su.neig_faces), deepcopy(su.neig_edges), deepcopy(su.face_area), deepcopy(su.edge_faces), deepcopy(su.face_normals))
#function used to equip a give TriSurface with its useful geometrical quantities
function genGeoQuant(surf::TriSurface)
	neig_faces=[[] for i in 1:surf.Nv]
	face_area=zeros(length(surf.faces))
	edges=[]
	face_normals=[[0.,0.,0.] for i in 1:length(surf.faces)]
	for f in 1:length(surf.faces)
		for vf in surf.faces[f]
			push!(neig_faces[Int(vf)], f)
		end
		ab=surf.vertices[surf.faces[f][1]].-surf.vertices[surf.faces[f][2]]
		ac=surf.vertices[surf.faces[f][3]].-surf.vertices[surf.faces[f][2]]
		face_area[f]=area(ab,ac)
		el=surf.faces[f]
		sg1=[el[1],el[2]]
		sg2=[el[2],el[3]]
		sg3=[el[1],el[3]]
		sg4=[el[2],el[1]]
		sg5=[el[3],el[2]]
		sg6=[el[3],el[1]]
		push!(edges, sg1)
		push!(edges, sg2)
		push!(edges, sg3)
		push!(edges, sg4)
		push!(edges, sg5)
		push!(edges, sg6)
		facen=cross_prod(ab,ac)
		facen/=norm(facen)
		face_normals[f]=facen
	end
	edges = filter(x->x[1]<x[2],edges)
	edges=union(edges)
	neig_edges=[[] for i in 1:surf.Nv]
	edge_faces=[[] for i in 1:length(edges)]
	for e in 1:length(edges)
		push!(neig_edges[edges[e][1]], e)			
		push!(neig_edges[edges[e][2]], e)
		for f in 1:length(collect(surf.faces[neig_faces[edges[e][1]][:]]))
			if edges[e][1] in collect(surf.faces[neig_faces[edges[e][1]][f]]) && edges[e][2] in collect(surf.faces[neig_faces[edges[e][1]][f]])
				push!(edge_faces[e], neig_faces[edges[e][1]][f])
			end
		end
	end
	
	return SurfGeoQuant(deepcopy(surf.Nv), deepcopy(surf.vertices), deepcopy(surf.faces), edges, deepcopy(surf.neig), neig_faces, neig_edges, face_area, edge_faces, face_normals)
end
#update surfgeoquant after vertex move
function update_geoquantvm(su::SurfGeoQuant, i, xn)
	#su2=SurfGeoQuant(su.Nv, su.vertices, su.faces, su.edges, su.neig, su.neig_faces, su.neig_edges, su.face_area, su.edge_faces, su.face_normals)
	su2=deepcopy(su)
	su2.vertices[i]=xn
	for f in su2.neig_faces[i]
		ab=su2.vertices[su2.faces[f][1]].-su2.vertices[su2.faces[f][2]]
		ac=su2.vertices[su2.faces[f][3]].-su2.vertices[su2.faces[f][2]]
		su2.face_area[f]=area(ab,ac)
		facen=cross_prod(ab,ac)
		facen/=norm(facen)
		su2.face_normals[f]=facen
	end
	return su2
end
#update surfgeoquant after link move
function update_geoquantlm(su::SurfGeoQuant, e)
	#suppose check for overlapping and maximum edge length constrains has already been done
	#su2=SurfGeoQuant(su.Nv, su.vertices, su.faces, su.edges, su.neig, su.neig_faces, su.neig_edges, su.face_area, su.edge_faces, su.face_normals)
	su2=deepcopy(su)
	t1,t2=su2.edge_faces[e][1],su2.edge_faces[e][2] #triangles (faces) sharing edge e
	e1,e2=su2.edges[e][1],su2.edges[e][2]
	ne1=filter(x -> (x!= e1 && x!= e2), su2.faces[t1])[1]
	ne2=filter(x -> (x!= e1 && x!= e2), su2.faces[t2])[1]
	#update edge
	su2.edges[e]=[ne1,ne2]
	#check for counterclockwise order and update faces
	if chechkcounterclock([su2.vertices[e1], su2.vertices[ne1], su2.vertices[ne2]])
		su2.faces[t1]=[e1, ne1, ne2]
	else
		su2.faces[t1]=[ne1, e1, ne2]
	end 	
	if chechkcounterclock([su2.vertices[e2], su2.vertices[ne2], su2.vertices[ne1]])
		su2.faces[t2]=[e2, ne2, ne1]
	else
		su2.faces[t2]=[e2, ne1, ne2]
	end
	#update neig
	filter!(x -> x!=e2, su2.neig[e1])
	filter!(x -> x!=e1, su2.neig[e2])
	push!(su2.neig[ne2], ne1)
	push!(su2.neig[ne1], ne2)
	#update neig faces
	filter!(x -> x!=t1 && x!=t2, su2.neig_faces[e1])
	filter!(x -> x!=t1 && x!=t2, su2.neig_faces[e2])
	push!(su2.neig_faces[ne1], t1)
	push!(su2.neig_faces[ne1], t2)
	su2.neig_faces[ne1]=union(su2.neig_faces[ne1])
	push!(su2.neig_faces[ne2], t1)
	push!(su2.neig_faces[ne2], t2)	
	su2.neig_faces[ne2]=union(su2.neig_faces[ne2])
	#update neig eddges
	filter!(x -> x!=e, su2.neig_edges[e1])
	filter!(x -> x!=e, su2.neig_edges[e2])
	push!(su2.neig_edges[ne1], e)
	push!(su2.neig_edges[ne2], e)
	#update face area and face normals
	for f in [t1, t2]
		ab=su2.vertices[su2.faces[f][1]].-su2.vertices[su2.faces[f][2]]
		ac=su2.vertices[su2.faces[f][3]].-su2.vertices[su2.faces[f][2]]
		su2.face_area[f]=area(ab,ac)
		facen=cross_prod(ab,ac)
		facen/=norm(facen)
		su2.face_normals[f]=facen
	end
	#update edge faces
	#edge between e1 and e2 remains unchanged, and the same holds for e1-ne1 and e2-ne2 
	edge_e1_ne2=0
	for ed in su2.neig_edges[e1]
		if su2.edges[ed][1]==ne2 || su2.edges[ed][2]==ne2
			edge_e1_ne2=ed
		end
	end
	edge_e2_ne1=0
	for ed in su2.neig_edges[e2]
		if su2.edges[ed][1]==ne1 || su2.edges[ed][2]==ne1
			edge_e2_ne1=ed
		end
	end
	filter!(x -> x!=t2, su2.edge_faces[edge_e1_ne2])
	push!(su2.edge_faces[edge_e1_ne2], t1)
	filter!(x -> x!=t1, su2.edge_faces[edge_e2_ne1])
	push!(su2.edge_faces[edge_e2_ne1], t2)
	return su2
end

#used just to prepare a Meshes.jl plot that can be plotted through makie
function prep4plot(vertices, faces)
	nv = Point3[(el[1],el[2],el[3]) for el in vertices]
	connec = connect.([(el[1],el[2],el[3]) for el in faces],Ngon)
	msh=SimpleMesh(nv, connec)
	return msh
end
