## You should be able to get stated with:

```
$ git clone git@github.com:Omar-Elrefaei/Debugging-GeneralizedSylvester.git
$ cd Debugging-GeneralizedSylvester
$ julia --project --threads 4
julia> using Revise; includet("investigation.jl")
julia> df = scan_seeds(func=gss_test, seeds=1:10)
```

This will run `gss_test` 10 time with problem parameters as defined at the top of the `gss_test` function (n=m=order=4).

## Caching & Performance
I found out that `cond` is the most expensive operation in `gss_test` particularly for large matrices.
So I setup memorization (ghost function `_cond`) such that the second time you do `scan_seeds(func=gss_test, seeds=1:10)` with the same params, it uses the cached results.

This came to be when I was running with `1:1000`, then decide that I want to add a new column to the dataframe and run again. 

This workflow with n=m=4 at order 4 if pretty nice. But still fairly painful at order 5.
For order 5, comment out the last dataframe columns (`cond`) and all the lines that mention `matEq|mateq` (gsylv from MatrixEquations takes about 8s at that scale)

This is why I ran a long n=m=4, order=5 scan over 5000 seeds and uploaded the resultant df as a csv.
you can read it with `DataFrame(CSV.File("order5.csv"))`


## Dataframe output
each `gss_test` invocation return an `OrderedDict` that represents a row in the `df` DataFrame; 
> If you need to add or remove columns from the dataframe, that should be done in the `dict = OrderedDict(...` line right away.

The rows then get merged together in the parent loop at `scan_seeds`.
Originally that loop is supposed to simply represent the following
```julia
for seed in seeds
	res = gss_test(seed)
	append!(df, res, cols=:union)
end
```
but things got complicated with the introduction of threads, atomics, progress-bar. And the try-catch for better `CTRL-C` cancellation.
Needless to say, I got a little carried away....

## Analysis
I've been doing rudementary analysis like this
- Find all cases where our gss iterative solution didn't satisfy the original problem (with scaled tolerance)
	- `filter(:gss_vs_orgPrb_rtol => isequal(false), df)` 
- Then take that and sort by the norm of the diffrences to more easily see any patterns
	- `filter(:gss_vs_orgPrb_rtol => isequal(false), df) |> x -> sort(x, ["norm(GSS-OrgPrb)"])` 
