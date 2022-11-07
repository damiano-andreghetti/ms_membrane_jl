#Dynamically Triangulated Monte Carlo in Julia 1.6.7
#Packages Meshes,MeshViz,Makie (for plotting), JLD for fileIO
#creation of triangulated sphere from icosahedral inspired by https://github.com/JanisErdmanis/LaplaceBIE.jl/blob/master/examples/sphere.jl
#TODO
#Listed-Link-Cell ?sistemare
#parallelize MC sweep
#prepare in SurfGeoQuant Nf1 and Nf2 (wrt to lower indexx vertex)
#using Meshes, MeshViz
#import GLMakie as Mke
include("sphere.jl")
include("MC_sweep.jl")
ico = create_ico()
sph=subNtimes(2,ico)
sph=normEdges(sph)
L=100
CELLS= genCells(sph,L,L,L, L) #radius of sphere circa =2^k where k is number of subdivision, so create cells a bit bigger
sgq = genGeoQuant(sph)
part=[false for i in 1:length(sgq.faces)]
#part[1]=true	
#sgq= Nsweeps(sgq, CELLS, 0.1, sqrt(3), 0.5, 500, part, 0.5,100)
sgq, mis=NsweepsPE(sgq,CELLS, 0.1, sqrt(3),0.5, 10000, 500, "pe_fixed_sub2")	
println(mis)
"""
mesh=prep4plot(sgq.vertices, sgq.faces)
fig = Mke.Figure(resolution = (800, 400))
viz(fig[1,1], mesh, showfacets=true)
fig
""" 