
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
	pass_fail_dict = Dict()
	# Threads.@threads 
	@showprogress for seed in seeds
		res = test_repetedly(func=func, seed=seed, n_tries=n_tries, verbose=false)
		pass_fail_dict[seed] = res
	end

	inconsistent_seeds = filter(p->ismissing(p.second), pass_fail_dict) |> sort 
	passed_seeds = filter(p->isequal(true)(p.second), pass_fail_dict) |> sort
	failed_seeds = filter(p->isequal(false)(p.second), pass_fail_dict) |> sort
	# printing block
	println("\nPassing seeds: $(length(passed_seeds))/$(length(seeds))")
	0<length(passed_seeds)<15 ? println("  =>	", [p.first for p in passed_seeds]) : nothing
	println("Failing seeds: $(length(failed_seeds))/$(length(seeds))")
	0<length(failed_seeds)<15 ? println("  =>	", [p.first for p in failed_seeds]) : nothing
	println("Inconsistent Seeds: $(length(inconsistent_seeds))/$(length(seeds))")
	0<length(inconsistent_seeds)<15 ? println("  =>	", [p.first for p in inconsistent_seeds]) : nothing
	if length(failed_seeds)+length(inconsistent_seeds) == 0 print("✨ Nice, no faliures or inconsistencies ✨") end
	if length(passed_seeds) == length(seeds) println("✨ All $(length(seeds)) seeds Passed ✨") end
	return passed_seeds
end	

function consistent_or_not(list::Vector{Bool})
	if isempty(list) return nothing
	elseif all(isequal(true), list) return true
	elseif all(isequal(false), list) return false
	else return missing end
end



######### Wrapper functions to test and trigger faliure #########

function gss_test(seed)
	if seed != false Random.seed!(seed) end
	n1 = 4
	n2 = 4
	order = 3


	# n1 = 2
	# n2 = 2
	# order = 3
	a = randn(n1,n1)
	b = randn(n1,n1)
	c = randn(n2,n2)
	ws = GeneralizedSylvesterWs(n1,n1,n2,order)
	d = randn(n1,n2^order)
	a_orig = copy(a)
	b_orig = copy(b)
	c_orig = copy(c)
	d_orig = copy(d)

	# @show "chkpt1"
	d = generalized_sylvester_solver!(a, b, c, copy(d_orig), order, ws)
	# @show "chkpt2"
	
	_a = kron_power(c_orig, order)
	# @show "chkpt2.5"
	lhs1 = a_orig*d + b_orig*d*_a
	# @show "chkpt3"
	
	rhs1 = d_orig

	res = lhs1 ≈ rhs1

	A = (kron(I(n2^order),a_orig) + kron(kron_power(c_orig', order),b_orig))
	b = vec(d_orig)
	sol2 = reshape(A\b, n1,n2^order)
	res2 = d ≈ sol2
	
	# if seed == 74 || seed == 75
	# if seed == 2 || seed == 3 || seed == 9 || seed == 62
	if res2== false
	# 	@info "result is $res"
	# 	@show seed
	# 	# @show a
	# 	# @show b
	# 	# @show c
	# 	# @show d_orig
	# 	# @show d
		@show norm(lhs1-rhs1)
	# 	@show norm(lhs1)
	# 	@show norm(rhs1)
	# 	@show cond(rhs1)
	# 	@show cond(lhs1)
	# 	println()
	end
	results = [res, res2]
	if !all(results)
		println("seed: $seed, res1 is $res, and res2 is $res2")
	end


	return res2


end

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