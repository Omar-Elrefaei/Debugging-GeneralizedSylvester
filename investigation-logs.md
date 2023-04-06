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





# Log
