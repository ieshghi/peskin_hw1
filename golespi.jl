module golespi

function ΔΘ(n::Int,ΔN::Int,N0::Int,Nb1::Int,Na1::Int)
	return ΔN+log(((N0-(n+1))/(Nb1+N0-(n+1))/(n/(Na1+n))))
end

function αp(n::Int,γ::Float64,ΔN::Int,N0::Int,Nb1::Int,Na1::Int)
		t = ΔΘ(n,ΔN,N0,Nb1,Na1)
		return γ*t/(exp(t)-1)
end

function αm(n::Int,γ::Float64,ΔN::Int,N0::Int,Nb1::Int,Na1::Int)
		t = ΔΘ(n-1,ΔN,N0,Nb1,Na1)
		return γ*t/(1-exp(-t))
end

function timestep(n::Int,γ::Float64,pars::Array{Int,1})
		tp = -log(rand())/αp(n,γ,pars[1],pars[2],pars[3],pars[4])		
		tm = -log(rand())/αm(n,γ,pars[1],pars[2],pars[3],pars[4])

		nout = n + (tp>tm)-(tp<tm)

		return nout,min(tp,tm)
end

function evolve(ni::Int,γ::Float64,pars::Array{Int,1},nsteps::Int)
		xs = zeros(nsteps+1)
		ts = zeros(nsteps+1)
		xs[1] = ni
		x = ni
		for i = 1:nsteps
				xnew,dt = timestep(x,γ,pars)
				xs[i+1] = xnew
				ts[i+1] = ts[i]+dt
				x = xnew
		end
		
		return xs,ts
end
end
