using StatsBase
#to evaluate autocorrelation time
#requires StatsBase
function at(vr::Vector{Float64})
	a=autocor(vr,collect(0:length(vr)-1))
	tau_int=0.5
	#display(plot(a))
	for t in 1:length(a)
		tau_int+=(a[t]*(1-(t/length(a))))
	end
	return tau_int
end
#block average a serie of data
#requires StatsBase
function block_ave(vr, sz)
	Nb=Int(floor(length(vr)/sz))
	val=[]
	for i in 1:Nb
		sb=1+(i-1)*sz
		eb=sb+sz-1
		push!(val, mean(vr[sb:eb]))
	end
	return val
end
#logarithmically ranged values
function logrange(x1, x2, n)
	return [10^y for y in range(log10(x1), log10(x2), length=n)]
end

