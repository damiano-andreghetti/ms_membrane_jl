using StatsBase,JLD, Plots
include("dami_utilities.jl")
function run_script()
	Ks=[10.0,30.0,100.0,300.0]
	k=Ks[1]
	println("k inclusione=",k)
	Bs=[]
	Ds=[]
	d=10
	println("distance d=",d)
	alphas::Vector{Float64}=[]
	betas::Vector{Float64}=[]	
	for cn in 10000:1000:6000000
		println(cn)
		h=load("K="*string(k)*"/config_"*string(cn)*".jld")["membrane"]
		p=load("K="*string(k)*"/config_"*string(cn)*".jld")["particles"]
		#display(wireframe(h, title=string(cn)))
		#misura 14+15
		#r_1=(50+d,50)     50 perchè la particella è a (50,50)
		r=50+d
		u_x=(h[r+1,50]-h[r-1,50])/2.0
		u_y=(h[r,50+1]-h[r,50-1])/2.0
		alpha=u_x^2-u_y^2
		#misura 16+17
		#r_2=(50,50+d) r_1=(50+d,50+d)
		u_x1=(h[r+1,r]-h[r-1,r])/2.0
		u_y2=(h[50,r+1]-h[50,r-1])/2.0
		beta=u_x1*u_y2
		push!(alphas,alpha)
		push!(betas,beta)
	end
	save("alpha_beta_d=10.jld", "alphas", alphas, "betas", betas)
	#alphas=load("alpha_beta_d=10.jld")["alphas"]
	#betas=load("alpha_beta_d=10.jld")["betas"]
	display(plot([alphas, betas]))
	alphas=block_ave(alphas, 100)
	betas=block_ave(betas, 100)
	al=mean(alphas)
	be=mean(betas)
	al_err=std(alphas)/sqrt(length(alphas)-1)
	be_err=std(betas)/sqrt(length(betas)-1)
	println("alpha=", al, "  error=", al_err)
	println("beta=", be, "  error=", be_err)
	k0=10.00
	B=-(4*be-al)*2*pi*pi*k0*k0*d*d
	D=-(2*al-4*be)*4*pi*pi*k0*k0*d*d
	B_err=2*sqrt(al_err^2+16*be_err^2)*pi*pi*k0*k0*d*d#((4*be_err)+al_err)*2*pi*pi*k0*k0*d*d
	D_err=8*sqrt(al_err^2+4*be_err^2)*pi*pi*k0*k0*d*d#((4*be_err)+2*al_err)*4*pi*pi*k0*k0*d*d
	println("B=",B," std=",B_err)
	println("D=",D," std=",D_err)		
	#display(plot([alphas, betas]))
end

#e=load("K=10.0/measures.jld")["energy"]
#display(plot(e))

run_script()