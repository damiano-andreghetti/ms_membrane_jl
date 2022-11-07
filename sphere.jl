using Meshes #requires Meshes package for perp4plot function
using JLD
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
#WORKS ONLY IF ORIGIN IS INSIDE THE MESH
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
	su2=su
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
	#suppose check for overlapping and maximum edge length constraints has already been done, also check for neighbour>=3 and no pyramidal structures
	#su2=SurfGeoQuant(su.Nv, su.vertices, su.faces, su.edges, su.neig, su.neig_faces, su.neig_edges, su.face_area, su.edge_faces, su.face_normals)
	su2=su
	t1,t2=su2.edge_faces[e][1],su2.edge_faces[e][2]
	x=su2.faces[t1]
	y=su2.faces[t2]
    x_new=0
    y_new=0            #x_new. y_new elements not contained in the edge, that will form new edge
    x_ind=0 			#indices of first vertex forming old edge in x and second vertex forming old edge in y
    y_ind=0
    u=deepcopy(su2.faces[t1])
	v=deepcopy(su2.faces[t2])
    both=[]
    for i in x
        if i in y
            push!(both, i)
		end
	end
    for i in 1:length(x)
        if x[i]== both[1]
            x_ind=i
		end
        if !(x[i] in both)
            x_new=x[i]
		end
	end
    for j in 1:length(y)
        if y[j]== both[2]
            y_ind=j
		end
        if !(y[j] in both)
            y_new=y[j]
		end
	end
    u[x_ind]=y_new
    v[y_ind]=x_new
	#update egdes (keeping them like (i,j) with i<j)
	if x_new<y_new
		su2.edges[e]=[x_new,y_new]
	else
		su2.edges[e]=[y_new,x_new]
	end
	#update faces
	su2.faces[t1]=u
	su2.faces[t2]=v
	#update neig
	b1,b2=both[1],both[2]
	filter!(x -> x!=b1, su2.neig[b2])
	filter!(x -> x!=b2, su2.neig[b1])
	push!(su2.neig[x_new], y_new)
	push!(su2.neig[y_new], x_new)
	#update neig faces
	filter!(x -> x!=t1 && x!=t2, su2.neig_faces[b1])
	filter!(x -> x!=t1 && x!=t2, su2.neig_faces[b2])
	for f in [t1,t2]
		if b1 in su2.faces[f]
			push!(su2.neig_faces[b1], f)
		end
		if b2 in su2.faces[f]
			push!(su2.neig_faces[b2], f)
		end
	end
	push!(su2.neig_faces[x_new], t1)
	push!(su2.neig_faces[x_new], t2)
	su2.neig_faces[x_new]=union(su2.neig_faces[x_new])
	push!(su2.neig_faces[y_new], t1)
	push!(su2.neig_faces[y_new], t2)	
	su2.neig_faces[y_new]=union(su2.neig_faces[y_new])
	#update neig_edges
	filter!(x -> x!=e, su2.neig_edges[b1])
	filter!(x -> x!=e, su2.neig_edges[b2])
	push!(su2.neig_edges[x_new], e)
	push!(su2.neig_edges[y_new], e)
	#update face area and normals
	for f in [t1, t2]
		ab=su2.vertices[su2.faces[f][1]].-su2.vertices[su2.faces[f][2]]
		ac=su2.vertices[su2.faces[f][3]].-su2.vertices[su2.faces[f][2]]
		su2.face_area[f]=area(ab,ac)
		facen=cross_prod(ab,ac)
		facen/=norm(facen)
		su2.face_normals[f]=facen
	end
	#update edge faces
	for vert in union([u; v]) #for all involved vertices
		for ed in su2.neig_edges[vert] #for all its neighbour edges
			su2.edge_faces[ed]=[]
			for f in 1:length(collect(su2.faces[su2.neig_faces[su2.edges[ed][1]][:]]))
				if su2.edges[ed][1] in collect(su2.faces[su2.neig_faces[su2.edges[ed][1]][f]]) && su2.edges[ed][2] in collect(su2.faces[su2.neig_faces[su2.edges[ed][1]][f]])
					push!(su2.edge_faces[ed], su2.neig_faces[su2.edges[ed][1]][f])
				end
			end
			su2.edge_faces[ed]=union(su2.edge_faces[ed])
		end
	end
	return su2
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
#same but for surfgeoquant
function genCells(surf::SurfGeoQuant, Nx,Ny,Nz, L)
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
#OTHER FUNCTIONS
#
#

#save configuration to file
function save_config(surf,filename)
	save("mesh_snapshot/"*filename*".jld", "surf", surf)
end
#load configuration from file
function load_config(filename)
	return load(filename)["surf"]
end

#used just to prepare a Meshes.jl plot that can be plotted through makie
function prep4plot(vertices, faces)
	nv = Point3[(el[1],el[2],el[3]) for el in vertices]
	connec = connect.([(el[1],el[2],el[3]) for el in faces],Ngon)
	msh=SimpleMesh(nv, connec)
	return msh
end
