using LinearAlgebra, Random

using KroneckerTools, QuasiTriangular, GeneralizedSylvesterSolver
using MatrixEquations, DataFrames, DataStructures

# Utilities for running long multi-threaded tests
using ProgressMeter, Memoization, ThreadSafeDicts

function scan_seeds(;func, seeds::UnitRange{Int64}, n_tries=1)
	lk = ReentrantLock()
	df = DataFrame()
	breakout = false

	prog = Progress(length(seeds))
	Threads.@threads for seed in seeds
		if breakout break end
		try
			res = gss_test(seed)
			lock(lk) do
				append!(df, res, cols=:union)
			end
			next!(prog)
		catch err 
			if err isa InterruptException breakout=true
			else rethrow() end
		end
	end

	if breakout @error "\nInterruptException but returning nevertheless" end
	passed_seeds = subset(df, :gss_vs_orgPrb) 
	failed_seeds = subset(df, :gss_vs_orgPrb=>x->x.==false)
	println("\nPassing seeds: $(length(eachrow(passed_seeds)))/$(length(seeds))")
	println("Failing seeds: $(length(eachrow(failed_seeds)))/$(length(seeds))")
	0<length(eachrow(failed_seeds))<15 ? println("  =>	", failed_seeds.seed) : nothing
	return df
end	

function gss_test(seed)
	if seed != false Random.seed!(seed) end
	n1 = 4
	n2 = 4
	order = 4
	rtol = sqrt(eps(1.0))*order
	rtol2 = sqrt(eps(1.0))*n1^order
	a = randn(Float64,n1,n1); b = randn(Float64,n1,n1); c = randn(Float64,n2,n2); d = randn(Float64,n1,n2^order)
	ws = GeneralizedSylvesterWs(n1,n1,n2,order)
	a_orig = copy(a); b_orig = copy(b); c_orig = copy(c); d_orig = copy(d)
	d = generalized_sylvester_solver!(a, b, c, copy(d_orig), order, ws)
	
	c_nKrons = kron_power(c_orig, order)
	lhs1 = b_orig*d*c_nKrons + a_orig * d 
	gss_vs_orgPrb = isapprox(lhs1, d_orig)
	gss_vs_orgPrb_rtol = isapprox(lhs1, d_orig, rtol=rtol)

	A = (kron(I(n2^order),a_orig) + kron(kron_power(c_orig', order),b_orig))
	sol2 = reshape(A\vec(d_orig), n1,n2^order)
	# sol2 = reshape( solve(LinearProblem(A, b)).u, n1,n2^order)
	res2 = isapprox(d, sol2)
	res3 = isapprox(d, sol2, rtol=rtol)
	res4 = isapprox(d, sol2, rtol=rtol2)
	# A*_s ≈ b ? nothing : @show A*_s ≈ b

	# sub-back the analytical solution into the original problem
	lhs2 = b_orig*sol2*c_nKrons + a_orig * sol2 
	res_ana_vs_orgPrb = isapprox(lhs2, d_orig)

	d_mateq = gsylv(b_orig, kron_power(c_orig,order), a_orig, 1.0, d_orig)
	res_mateq_vs_gss = isapprox(d, d_mateq)
	res_mateq_vs_ana = isapprox(sol2, d_mateq)
	
	lhs3 = b_orig*d_mateq*c_nKrons + a_orig * d_mateq 
	res_matEq_vs_orgPrb = isapprox(lhs3, d_orig)
	dict = OrderedDict(
		:seed => seed,
		:gss_vs_orgPrb => gss_vs_orgPrb,
		:gss_vs_orgPrb_rtol => gss_vs_orgPrb_rtol,
		:gss_vs_ana => res2,
		:gss_vs_ana_rtol => res3,
		:gss_vs_ana_rtol2 => res4,
		:ana_vs_orgPrb => res_ana_vs_orgPrb,
		:matEq_vs_gss => res_mateq_vs_gss,
		:matEq_vs_ana => res_mateq_vs_ana,
		:matEq_vs_orgPrb => res_matEq_vs_orgPrb,
		Symbol("norm(GSS-OrgPrb)") => norm(lhs1-d_orig),
		Symbol("norm(GSS-Ana)") => norm(d-sol2),
		Symbol("norm(Ana-OrigProb)") => norm(lhs2-d_orig),
		Symbol("norm(matEq-GSS)") => norm(d-d_mateq),
		Symbol("norm(matEq-Ana)") => norm(sol2-d_mateq),
		Symbol("norm(matEq-OrigProb)") => norm(lhs3-d_orig),
		Symbol("cond(Ana)") => _cond(A),
		Symbol("cond(b&d)") => _cond(b)+_cond(d_orig),
		)

	return dict
end

@memoize ThreadSafeDict function _cond(x)
	return cond(x) 
end


######### util from runtest
function kron_power(x,order)
    if order == 0
        return 1
    elseif order == 1
        return x
    else
        m, n = size(x)
        y = Matrix{Float64}(undef, m^order, n^order)
        y1 = similar(y)
        v = view(y, 1:m, 1:n)
        v1 = view(y1, 1:m, 1:n)
        v .= x
        for i = 1:(order-1)
            tmp = v1
            v1 = v
            v = tmp
            v = view(v.parent, 1:m^(i + 1), 1:n^(i + 1))
            kron!(v, v1, x)
        end
    end
    return v.parent
end