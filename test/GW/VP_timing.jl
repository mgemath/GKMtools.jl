n_repeats = 5

F = gkm_3d_twisted_flag();
beta = curve_class(F, "3", "4");

d = 5
g = 4
@btime gw = gromov_witten(F, d * beta, 0, class_one(); g=g, show_bar=false)

# usando _2
# # 448.284 s (474314692 allocations: 23.77 GiB)
# 1125103//14400

# usando _3 
#   419.109 s (474313283 allocations: 23.76 GiB)
# 1125103//14400

@btime begin
  for d in 1:4
    for g in 1:3
      gw = gromov_witten(F, d * beta, 0, class_one(); g=g, show_bar=false)
      println("Genus $g, degree $d: $gw")
    end
  end
end

# We run the above code 3-4 times for the above code without psi classes.
# the goal is to see if it makes sense to store vertex polynomials or if H is enough.

######

# n_repeats = 5, d in 1:3, g in 1:3, twisted flag class_one().

# with VPs:
# 	59.962779 seconds (18.95 M allocations: 918.598 MiB, 0.72% gc time, 7.78% compilation time)
# 	59.121890 seconds (28.04 M allocations: 1.221 GiB, 1.13% gc time, 6.90% compilation time)
# 	60.185003 seconds (28.04 M allocations: 1.221 GiB, 1.13% gc time, 6.69% compilation time)

# with H:
# 	61.457161 seconds (30.84 M allocations: 1.357 GiB, 1.20% gc time, 9.58% compilation time: 18% of which was recompilation)
# 	62.782081 seconds (28.04 M allocations: 1.222 GiB, 1.01% gc time, 6.58% compilation time)
# 	59.600786 seconds (28.04 M allocations: 1.223 GiB, 1.09% gc time, 6.63% compilation time)
# 	60.591393 seconds (28.04 M allocations: 1.223 GiB, 1.08% gc time, 6.56% compilation time)

#######

# repeats = 5, d in 1:4, g in 1:3, twisted flag class_one().

# With VPs:
# 	399.041496 seconds (156.94 M allocations: 6.818 GiB, 0.92% gc time, 1.54% compilation time)

# With H:
# 	377.038860 seconds (157.99 M allocations: 6.866 GiB, 0.96% gc time, 1.79% compilation time)