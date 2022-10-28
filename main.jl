#Dynamically Triangulated Monte Carlo in Julia 1.6.7
#Packages Meshes,MeshViz,Makie (for plotting)
#creation of triangulated sphere from icosahedral taken from by https://github.com/JanisErdmanis/LaplaceBIE.jl/blob/master/examples/sphere.jl
#TODO
#Listed-Link-Cell ?sistemare
#parallelize MC sweep
#prepare in SurfGeoQuant Nf1 and Nf2 (wrt to lower indexx vertex)
using Meshes, MeshViz
import GLMakie as Mke
include("sphere.jl")
include("MC_sweep.jl")
ico = create_ico()
sph=subNtimes(3,ico)
sph=normEdges(sph)
CELLS= genCells(sph,3,3,3,25) #radius of sphere circa =2^k where k is number of subdivision, so create cells a bit bigger
sgq = genGeoQuant(sph)
part=[false for i in 1:length(sgq.faces)]
part[1]=true	
#@time ev_energy(sgq, 500, part, 0.5)
sgq = Nsweeps(sgq, CELLS, 0.1, sqrt(3), 0.5, 5, part, 3,8, 1000)	

mesh=prep4plot(sgq.vertices, sgq.faces)
fig = Mke.Figure(resolution = (800, 400))
viz(fig[1,1], mesh, showfacets=true)
fig
