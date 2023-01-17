#Dynamically Triangulated Monte Carlo in Julia 1.8.2
#Packages Meshes,MeshViz,Makie (for plotting), JLD for fileIO
#creation of triangulated sphere from icosahedral inspired by https://github.com/JanisErdmanis/LaplaceBIE.jl/blob/master/examples/sphere.jl
using JLD
using Plots
using Statistics
using Meshes, MeshViz
import GLMakie as Mke
gr()
include("sphere.jl")
include("MC_sweep.jl")
println("libraries loaded")
ico = create_ico()
sph=subNtimes(3,ico)
sph=normEdges(sph)
L=100
CELLS= genCells(sph,L,L,L, L) #radius of sphere circa =2^k where k is number of subdivision, so create cells a bit bigger
sgq = genGeoQuant(sph)
println(norm(sgq.vertices[1]))
PartIns=64
part=Int16[-1 for i in 1:sgq.Nv]
#part[shuffle(collect(1:sgq.Nv))[1:PartIns]].=+1
#part[1:12].=1
sigma=0.07
lmax=sqrt(3)
radius=0.5
mu=0.0
k=6.0
Nsvert=sgq.Nv
Nslink=3*sgq.Nv
Nspart=0#PartIns
Nsweep=50000
#kvals=[1,5,10,20,50,100,200,500]
#sigmavals=[0.3,0.07,0.04,0.023,0.017,0.01,0.006, 0.003]
#these are to have approx 50% acceptance
function sweeprun(sgq,CELLS,sigma, lmax, radius, k,part, mu, Nsvert,Nslink,Nspart,N)
	x=[]
	y_rad=[]
	y_acc=[]
	for l in 1:N
		if l%100==0
			println("step ", l)
		end
		if l%1000==0
			save_config(sgq,part,"k="*string(k)*"/config_Nv642_k="*string(k)*"_step"*string(l))
		end
		sgq, part, a=MC_sweepPM(sgq,CELLS,sigma, lmax, radius, k,part, mu,Nsvert,Nslink,Nspart, verbose=false)
		"""
		if l>14000 && l <16000 && l%100==0
			mesh=prep4plot(sgq.vertices, sgq.faces)
			fig = Mke.Figure(resolution = (800, 400))
			viz(fig[1,1], mesh, showfacets=true)
			for i in 1:sgq.Nv
				magn=0.3
				if part[i]==+1
					viz!(fig[1,1],Sphere((sgq.vertices[i][1],sgq.vertices[i][2],sgq.vertices[i][3]),radius*magn), color=:red)
				end
			end
			Mke.save("frames/step_"*string(l)*".png",fig)
			fig=Nothing
		end
		"""
		push!(x,l)
		if l>20000
			push!(y_rad, gyr_rad(sgq))
		end
		push!(y_acc, a[1])
	end
	return x,y_rad, y_acc
end
#sgq,part=load_config("k="*string(k)*"/equipartition_Nv162_k="*string(k)*"_step1000000.jld")
#part=convert(Vector{Int16}, part)
#CELLS=genCells(sgq,L,L,L,L)
a,b,c=sweeprun(deepcopy(sgq),deepcopy(CELLS),sigma, lmax, radius, k,deepcopy(part), mu, Nsvert,Nslink,Nspart,Nsweep)
b=convert(Vector{Float64}, b)
save("k="*string(k)*"/gyr_rad_square.jld","Rg",b)

