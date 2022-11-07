#plot config
using Meshes, MeshViz
import GLMakie as Mke
include("sphere.jl")
sgq=load_config("mesh_snapshot/pe_fixed_sub2_step3500.jld")
println("mesh loaded")
mesh=prep4plot(sgq.vertices, sgq.faces)
println("prepared for plot")
fig = Mke.Figure(resolution = (800, 400))
viz(fig[1,1], mesh, showfacets=true)
println("plotting")
fig
