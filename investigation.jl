
function test_repetedly(;func, seed=false, n_tries=20, verbose=true)
	pass_fail_log::Vector{Bool} = []
	for i in 1:n_tries
		try
			r = func(seed)
			push!(pass_fail_log, r)		
		catch ex
			if ex == SingularException
				if seed != false return nothing # if seeded, skip the rest
				else continue end
			else rethrow()
		end
	end
	end
	# true, false, or missing
	results = consistent_or_not(pass_fail_log) 
	if verbose
		if ismissing(results) println("Damn! Some failed but some passed ❕")
		elseif results println("Passed ✔") 
		elseif !results println("Failed ❌")
		end
	end
	return results
end

function scan_seeds(;func, seeds::UnitRange{Int64}, n_tries=1)
	lk = ReentrantLock()
	pass_fail_dict = Dict()
	df = DataFrame()

	prog = Progress(length(seeds))
	Threads.@threads for seed in seeds
		res = gss_test(seed)
		lock(lk) do
			append!(df, res, cols=:union)
		end
		next!(prog)
	end

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
	
	_a = kron_power(c_orig, order)
	lhs1 = b_orig*d*_a + a_orig * d 
	gss_vs_orgPrb = isapprox(lhs1, d_orig)
	gss_vs_orgPrb_rtol = isapprox(lhs1, d_orig, rtol=rtol)

	A = (kron(I(n2^order),a_orig) + kron(kron_power(c_orig', order),b_orig))
	sol2 = reshape(A\vec(d_orig), n1,n2^order)
	# sol2 = reshape( solve(LinearProblem(A, b)).u, n1,n2^order)
	res2 = isapprox(d, sol2)
	res3 = isapprox(d, sol2, rtol=rtol)
	res4 = isapprox(d, sol2, rtol=rtol2)
	# A*_s ≈ b ? nothing : @show A*_s ≈ b

	## subbing in the linAlg solution into the original problem
	lhs2 = b_orig*sol2*_a + a_orig * sol2 
	res_ana_vs_orgPrb = isapprox(lhs2, d_orig)

	if res2 == false
		# println("seed: $seed, res1 is $res, and res2 is $res2")
		# @show lhs1, d_orig, lhs2, d, sol2
		# if abs(cond(A)) > 1e6
		# 	@show "false: cond number is $(cond(A))"
		# end
		# @printf "seed is %s:  " seed
		@printf "norm(lhs1-d_orig) = %.1E %s"        norm(lhs1-d_orig) res
		@printf "\t-  norm(d-sol2) = %.1E %s"        norm(d-sol2) res2
		@printf "\t-  norm(lhs2-d_orig) = %.1E %s\n" norm(lhs2-d_orig) res3
	end

	d_mateq = gsylv(b_orig, kron_power(c_orig,order), a_orig, 1.0, d_orig)
	res_mateq_vs_gss = isapprox(d, d_mateq)
	res_mateq_vs_ana = isapprox(sol2, d_mateq)
	
	lhs3 = b_orig*d_mateq*_a + a_orig * d_mateq 
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
		Symbol("cond(Ana)") => cond(A),
		Symbol("cond(b&d)") => cond(b)+cond(d_orig),
		)

	# res3 is worst, and surprisingly not corolated with res2
	return dict

end

function consistent_or_not(list::Vector{Bool})
	if isempty(list) return nothing
	elseif all(isequal(true), list) return true
	elseif all(isequal(false), list) return false
	else return missing end
end

######### Wrapper functions to test and trigger faliure #########


function linsolv_test(seed)
	if seed != false Random.seed!(seed) end
	n1 = 4
	n2 = 3
	order = 2
	a = randn(n1,n1)
	b = randn(n1,n1)
	c = randn(n2,n2)
	ws = GeneralizedSylvesterWs(n1,n1,n2,order)
	d = randn(n1,n2^order)
	a_orig = copy(a)
	b_orig = copy(b)
	c_orig = copy(c)
	d_orig = copy(d)

	d = generalized_sylvester_solver!(a, b, c, copy(d_orig), order, deepcopy(ws))
	a_orig*d + b_orig*d*kron(c_orig,c_orig) ≈ d_orig
end



# extracted from runtests and scanned to find faliure
function solve1_test(seed)
	if seed != false Random.seed!(seed) end
	n = 4
	m = 4
	order = 3
    # __ = randn(m*n^order)
	a = randn(m,m)  ## matters here the first one at 3,5,5,1
    b = randn(m,m)
    c = randn(n,n)
    d = randn(m*n^order)
    alpha = randn()
    beta = randn()

    #only real eigenvalues
    # t = lu(a\b).U
   	# s = lu(c).U
   
    t = schur(a\b).T
    s = schur(c).T

    t2 = QuasiUpperTriangular(t*t)
    s2 = QuasiUpperTriangular(s*s)
    t = QuasiUpperTriangular(t)
    s = QuasiUpperTriangular(s)
    r = 1.0
    nd = m*n^(order-1)
    
    ws = GeneralizedSylvesterWs(m, m, n, order)
	d = randn(m*n^order)
	d_orig = copy(d)

	GeneralizedSylvesterSolver.solve1!(r,order,t,t2,s,s2,d,ws)
	
	lhs = (I(m*n^order) + kron(kron_power(s',order),r*t))
	lhs_sp = SparseMatrixCSC(lhs)
	d_target = solve(LinearProblem(lhs_sp, d_orig)).u
    # d_target = lhs \ d_orig

	res = d ≈ d_target
	if res == false
		# cond_num = abs(cond(lhs)) # 9 is too lax, could go 10 (np 11 is too tight sometimes lets 1/1000)
		# cond_num = abs(condskeel(lhs)) # 1e9 (1e8 to catch seed 1288 but let way too many false-pos)		
		cond_num = abs(cond(lhs)) 
		if cond_num > 1e8
	        @info "Seed $seed: condition number $cond_num is too high, skipping..."
	        throw(SingularException)
		else
			@warn "Seed $seed: is false but condition number is $cond_num"
		end
    end

    return res
end


function solveiip_test(seed)
	if seed != false Random.seed!(seed) end
	order = 2
	m=6
	n=6
    a = randn(m,m)
    b = randn(m,m)
    c = randn(n,n)
    d = randn(m*n^order)

    alpha = randn()
    beta = randn()
    t = schur(a\b).T
    s = schur(c).T
    t2 = QuasiUpperTriangular(t*t)
    s2 = QuasiUpperTriangular(s*s)
    t = QuasiUpperTriangular(t)
    s = QuasiUpperTriangular(s)
    r = 1.0
    nd = m*n^(order-1)
     ws = GeneralizedSylvesterWs(m, m, n, order)

 
    d = randn(2*m*n^(order-1))
    d_orig = copy(d)
    a = randn()
    b1 = -abs(randn())
    b2 = abs(randn())
    G = [a b1; b2 a]
    GeneralizedSylvesterSolver.solviip2(alpha,beta,a,b2,b1,order,t,t2,s,s2,d,ws)
    lhs = (I(2*m*n^(order-1)) + 2*alpha*kron(kron(G',kron_power(s',order-1)),t)
                + (alpha*alpha + beta*beta)*kron(kron(G'*G',kron_power(s2',order-1)),t2))
    
    # @show cond(lhs)

    d_target = lhs\d_orig
    @show maximum(abs.(d .- d_target))
    d ≈ d_target

end


function pure_generalized_sylvester_test(seed)
	if seed != false Random.seed!(seed) end
	n1 = 4
	n2 = 4
	order = 2
	a = randn(n1,n1)
	b = randn(n1,n1)
	c = randn(n2,n2)
	d = randn(n1, n2^order)
	a_orig = copy(a)
	b_orig = copy(b)
	c_orig = copy(c)
	d_orig = copy(d)
	ws = GeneralizedSylvesterWs(n1,n1,n2,order)

	# AX + BXC = D
	# solving a x + b x (c ⊗ c ⊗ ... ⊗ c) = d
	# d = generalized_sylvester_solver!(a, b, c, copy(d_orig), order, deepcopy(ws))
	

	# X = gsylv(A,B,D,y,E)....solves AXB + DXy = E
	# d = gsylv(b, c, a, 1.0, d_orig)
	
	d = gsylv(b, kron(c,c), a, 1.0, d_orig)


	lhs1 = a_orig*d + b_orig*d*kron(c_orig,c_orig)
	rhs1 = d_orig
	lhs1 ≈ rhs1 

	# d ≈ d0
end


function addition_demo_test(seed)
	if seed != false Random.seed!(seed) end
	a = randn()
	b = randn()
	+(a,b) ≈ a+b
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