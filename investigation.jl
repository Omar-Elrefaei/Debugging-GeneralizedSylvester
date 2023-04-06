function test_repetedly(;func, seed=false, n_tries=20, verbose=true)
	pass_fail_log::Vector{Bool} = []
	for i in 1:n_tries
		push!(pass_fail_log, func(seed))
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

function scan_seeds(;func, seeds::UnitRange{Int64}, n_tries=20)
	pass_fail_dict = Dict()
	Threads.@threads for seed in seeds
		res = test_repetedly(func=func, seed=seed, n_tries=n_tries, verbose=false)
		pass_fail_dict[seed] = res
	end

	# print the seeds that consistently passed or was inconsistent, the rest has failed
	inconsistent_seeds = filter(p->ismissing(p.second), pass_fail_dict) |> sort 
	passed_seeds = filter(p->isequal(true)(p.second), pass_fail_dict) |> sort
	failed_seeds = filter(p->isequal(false)(p.second), pass_fail_dict) |> sort
	println("Passed seeds: ", [p.first for p in passed_seeds])
	println("Failed seeds: ", [p.first for p in failed_seeds])
	println("Inconsistent seeds: ", [p.first for p in inconsistent_seeds])
	if length(passed_seeds) == length(seeds) print("✨ All $(length(seeds)) seeds Passed ✨") end
	return passed_seeds
end


function consistent_or_not(list::Vector{Bool})
	if all(isequal(true), list) return true
	elseif all(isequal(false), list) return false
	else return missing end
end


function gss_test(seed)
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
	order = 2
	n = 6
	m = 6
    # __ = randn(m*n^order)
	a = randn(m,m)  ## matters here the first one at 3,5,5,1
    b = randn(m,m)
    c = randn(n,n)
    d = randn(m*n^order)
    alpha = randn()
    beta = randn()
    t = lu(a\b).U
    s = lu(c).U
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
	d_target = (I(m*n^order) + kron(kron_power(s',order),t))\d_orig
	d ≈ d_target
    
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