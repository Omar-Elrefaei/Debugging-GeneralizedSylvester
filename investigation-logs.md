- confirmed ldiv does the right thing
- both Kron calls call code that is tested at 2nd order

# Weridest bug I've ever faced
- going back to old deps, GSS is passing...no failing. BAH intermittent
- consistent pass/fail depending on Random.seed()
- Construst a testing setup before digging into solv*
...
- 2H later now we have a nice setup
- was able to trigger faliures in solve1
	- Little hard to trigger at orders lower than 3, but found at order=2, n=4, m=3 (seeds 652 and 727)
	- order 2 with n, m 6 or up is easy
...
...
*Heavy call with Michel*
*After lots of debugging, he decided to investigate cond(lhs)*
- Turns out yes, lhs in solve1 "test" is sometimes close to singular so inverting it is not very accurate
**Conclusion:** we assume the problem is with the testing ldiv a not-very-inverable matrix and not our functions 
- but after seemingly solving the problem 2 things comeup (one: not very easy to filter out all true-poisives without a ton of false-positives...this doesn't seem to be exactly the right createria but close)
- Two: that doesn't really solve the problem with generalized_sylvester_solver! since there is no ldiv there
	(with solve1 the \approx values where actually close, that is not the case with gss)
	- I have few ideas: 1: intercept args to solve1 call in a failing gss and see if solve1 is actually failing there
	- replace the solve1 call in gss with a julia ldiv version like the tests
	- 
- Observation: gss, number of failing seeds increase with n (equality faliure is big)





replacing norm with opnorm in isapprox
pdiv instead of ldiv



cond2 -> opnorm2 -> svdvals
> Note in pinv:  For inverting dense ill-conditioned matrices in a least-squares sense, rtol =
  sqrt(eps(real(float(one(eltype(M)))))) is recommended.

# Log
