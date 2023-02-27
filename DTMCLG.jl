#
#
#Here lattice gas dynamics is implemented on triangulated surface (defined in sphere.jl), which evolves through 
#Metropolis Monte Carlo (defined in MC_sweep.jl) 
#
#
include("sphere.jl")
include("MC_sweep.jl")
function part_sweep(sgq::SurfGeoQuant,part::Vector{Bool},Ki::Float64,Kd::Float64,g::Float64,Ne::Int)
	for v in shuffle(collect(1:sgq.Nv))
		if part[v]==0 && rand() < Ki
			#site empty, try to inject particle
			part[v]=1
		elseif part[v]==1 
			empty_neig=Int[]
			for j in sgq.neig[v]
				if part[j]==0
					push!(empty_neig, j)
				end
			end
			occupied::Int=length(sgq.neig[v])-length(empty_neig)
			if length(empty_neig)>0 && rand()<(Kd/(g^occupied))
				part[v]=0
				part[empty_neig[rand(1:length(empty_neig))]]=true
			end
		end
	end
	clusters=find_clusters(sgq,part)
	#extract clusters larger than Ne
	for cl in clusters
		if length(cl)>=Ne
			for p in cl
				part[p]=0
			end
		end
	end
	return part
end
function runDTMCLG(sgq::SurfGeoQuant,CELLS::CellsGrid,l0::Float64, lmax::Float64, radius::Float64, k::Float64,sigma::Float64,part::Vector{Bool}, A::Float64,B::Float64, Nsvert::Int,Nslink::Int,TMC::Int,Nrun::Int,Ninf::Int,Ki::Float64,Kd::Float64,g::Float64,Ne::Int, folder::String)
	x=[]
	y=[]
	Ene=ev_energy(sgq,k,sigma,part,A,B) #here will be stored evaluated energies, to avoid useless computations
	#save parameters in the folder for the simulation
	save(folder*"/parameters.jld", "l0",l0,"lmax",lmax,"radius",radius,"k",k,"sigma",sigma,"A",A,"B",B,"Nsvert",Nsvert,"Nslink",Nslink,"TMC",TMC,"Ki",Ki,"Kd",Kd,"g",g,"Ne",Ne)
	for l in 1:Nrun	
		if l%Ninf==0
			println("step ", l)
			save_config_fastIO(sgq,part,folder*"/config_step"*string(l))
			"""
			mesh=prep4plot(sgq.vertices, sgq.faces)
			fig = Mke.Figure(resolution = (800, 400))
			viz(fig[1,1], mesh, showfacets=true)
			for i in 1:sgq.Nv
				magn=0.4
				if part[i]==+1
					viz!(fig[1,1],Sphere((sgq.vertices[i][1],sgq.vertices[i][2],sgq.vertices[i][3]),radius*magn), color=:red)
				end
			end
			Mke.save("frames/step_"*string(l)*".png",fig)
			fig=Nothing
			"""
		end
		for MCs in 1:TMC
			sgq, part,Ene=MC_sweepPM(sgq,CELLS,Ene,l0, lmax, radius, k,sigma,part, A,B,Nsvert,Nslink, verbose=false)
		end
		part=part_sweep(sgq,part,Ki,Kd,g,Ne)
		push!(x,l)
		push!(y, sum(Ene))
	end
	return x,y
end

