module golespi
using Plots; pyplot()
using Statistics

function ΔΘ(n::Int,ΔN::Int,N0::Int,Nb1::Int,Na1::Int)
	return ΔN+log(((N0-(n+1))/(Nb1+N0-(n+1))/(n/(Na1+n))))
end

function αp(n::Int,ΔN::Int,N0::Int,Nb1::Int,Na1::Int)
		t = ΔΘ(n,ΔN,N0,Nb1,Na1)
		return t/(exp(t)-1)
end

function αm(n::Int,ΔN::Int,N0::Int,Nb1::Int,Na1::Int)
		t = ΔΘ(n-1,ΔN,N0,Nb1,Na1)
		return t/(1-exp(-t))
end

function timestep(n::Int,pars::Array{Int,1})
		tp = -log(rand())/αp(n,pars[1],pars[2],pars[3],pars[4])		
		tm = -log(rand())/αm(n,pars[1],pars[2],pars[3],pars[4])

		nout = n + (tp>tm)-(tp<tm)

		return nout,min(tp,tm)
end

function eqval(pars::Array{Int,1})
		ΔN = pars[1]
		N0 = pars[2]
		Nb1 = pars[3]
		Na1 = pars[4]

		if ΔN>0
			return (1-Nb1-N0+exp(ΔN)*(N0-Na1-1)+sqrt(4*Na1*exp(ΔN)*(exp(ΔN)-1)*(N0-1)+(Nb1-1+exp(ΔN)*(1+Na1-N0)+N0)^2))/(2*(exp(ΔN)-1))
		else 
			return Na1*(N0-1)/(Na1+Nb1)
		end
end

function evolve(ni::Int,pars::Array{Int,1},nsteps::Int)
		xs = zeros(nsteps+1)
		ts = zeros(nsteps+1)
		xs[1] = ni
		x = ni
		for i = 1:nsteps
				xnew,dt = timestep(x,pars)
				xs[i+1] = xnew
				ts[i+1] = ts[i]+dt
				x = xnew
		end
		eq = eqval(pars)	
		return xs,ts,eq
end

function hw_script()
	pars = [0,1000,700,300]
	nt = 10^5
	x,t,e = evolve(100,pars,nt)

	pars = [2,1000,700,300]
	nt = 10^5
	x2,t2,e2 = evolve(100,pars,nt)

	pars = [0,1000,600,400]
	nt = 10^5
	x3,t3,e3 = evolve(100,pars,nt)

	plot([t,t,t2,t2,t3,t3],[x,e*ones(size(x)),x2,e2*ones(size(x2)),x3,e3*ones(size(x3))],legend=false,reuse = false)
	
	nt = 10^6
	pa = 0.6
	pb = 0.4
	nvals = Int.(floor.(exp.(2.4*LinRange(2,4,20))))
	st = zeros(size(nvals))
	for i = 1:length(nvals)
		na = 10*Int(floor(pa*nvals[i]/10))
		nb = 10*Int(floor(pb*nvals[i]/10))
		nvals[i] = 10*Int(floor(nvals[i]/10))
		pars = [0,nvals[i],na,nb]
		x,t,e = evolve(100,pars,nt)
		st[i] = std(x[10^5:end])
	end
	
	plot([log.(nvals),log.(nvals[10:15])],[log.(nvals[10:15])/2,log.(st)])

	return nvals,st
end

end
