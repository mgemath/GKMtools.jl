# The point of this file is to look for smooth Schubert varieties that are not positive.

function search_nonpositive_Schubert(R::RootSystem; printFully::Bool=false)
	res = Vector{Any}()
	S = simple_roots(R)
	for s in 0:length(S)

		for S_sub in combinations(S, s)

			S_sub_vect = collect(S_sub)
			BO = get_bruhat_order_of_generalized_flag(R, S_sub_vect)

			for l in BO.labels
				Sch = generalized_gkm_schubert(R, S_sub_vect, l).self
				# Only consider rationally smooth Schubert varieties
				!isvalid(Sch; printDiagnostics=false) && continue
				# Only consider positive dimensional Schubet varieties
				isempty(edges(Sch.g)) && continue
				#println((R, S_sub_vect, l))

				min_chern_number = 0
				try
					min_chern_number = minimum(e -> chern_number(e, Sch), edges(Sch.g))
				catch e
					if isa(e, ArgumentError) # This happens if the Schubert variety is rationally smooth (i.e. all vertices have the same valency) but not smooth.
						#println("Argument error for $((R, S_sub_vect, l))")
						continue
					else
						rethrow()
					end
				end
				if min_chern_number <= 0
					push!(res, (R, S_sub_vect, l, min_chern_number))
					if printFully
						println()
						println("$((R, S_sub_vect, l)):")
						println("Schubert variety is: $Sch")
						print_curve_classes(Sch; printConAsForZeroC1=true)
					end
				end
			end
		end
	end
	return res
end

# julia> search_nonpositive_Schubert(root_system(:A, 3); printFully=true)

# (Root system of type A3, RootSpaceElem[], "s1*s3*s2*s1"):
# Schubert variety is: GKM graph with 12 nodes and valency 4
# s1 -> id: (0, 0, 1), Chern number: 2
# s2*s1 -> s1: (0, 1, 1), Chern number: 4
# s1*s2*s1 -> id: (0, 1, 1), Chern number: 4
# s1*s2*s1 -> s2*s1: (0, 1, 0), Chern number: 2
# s3*s2*s1 -> s2*s1: (1, 1, 1), Chern number: 4
# s1*s3*s2*s1 -> s1*s2*s1: (1, 1, 1), Chern number: 4
# s1*s3*s2*s1 -> s3*s2*s1: (0, 1, 0), Chern number: 2
# s2 -> id: (0, 1, 0), Chern number: 2
# s2 -> s2*s1: (0, 0, 1), Chern number: 2
# s1*s2 -> s1: (0, 1, 0), Chern number: 2
# s1*s2 -> s1*s2*s1: (0, 0, 1), Chern number: 2
# s1*s2 -> s2: (0, 1, 1), Chern number: 4
# s3*s2 -> s3*s2*s1: (0, 0, 1), Chern number: 2
# s3*s2 -> s2: (1, 1, 0), Chern number: 2
# s1*s3*s2 -> s1*s3*s2*s1: (0, 0, 1), Chern number: 2
# s1*s3*s2 -> s1*s2: (1, 1, 0), Chern number: 2
# s1*s3*s2 -> s3*s2: (0, 1, 1), Chern number: 4
# s3 -> id: (1, 0, 0), Chern number: 0, con-as: ZZRingElem[2, -1, -1, 0]
# s3 -> s1*s3*s2*s1: (0, 1, 1), Chern number: 4
# s3 -> s3*s2: (0, 1, 0), Chern number: 2
# s1*s3 -> s1: (1, 0, 0), Chern number: 0, con-as: ZZRingElem[2, -1, -1, 0]
# s1*s3 -> s3*s2*s1: (0, 1, 1), Chern number: 4
# s1*s3 -> s1*s3*s2: (0, 1, 0), Chern number: 2
# s1*s3 -> s3: (0, 0, 1), Chern number: 2

# (Root system of type A3, RootSpaceElem[], "s1*s2*s3*s2"):
# Schubert variety is: GKM graph with 12 nodes and valency 4
# s1 -> id: (0, -1, 1), Chern number: 0, con-as: ZZRingElem[2, -1, -1, 0]
# s2 -> id: (1, 0, 0), Chern number: 2
# s1*s2 -> s1: (1, 0, 0), Chern number: 2
# s1*s2 -> s2: (1, -1, 1), Chern number: 2
# s3*s2 -> s2: (0, 1, 0), Chern number: 4
# s1*s3*s2 -> s1*s2: (0, 1, 0), Chern number: 4
# s1*s3*s2 -> s3*s2: (1, -1, 1), Chern number: 2
# s2*s3*s2 -> id: (0, 1, 0), Chern number: 4
# s2*s3*s2 -> s3*s2: (-1, 1, 0), Chern number: 2
# s1*s2*s3*s2 -> s1: (0, 1, 0), Chern number: 4
# s1*s2*s3*s2 -> s1*s3*s2: (-1, 1, 0), Chern number: 2
# s1*s2*s3*s2 -> s2*s3*s2: (0, 0, 1), Chern number: 4
# s3 -> id: (-1, 1, 0), Chern number: 2
# s3 -> s3*s2: (1, 0, 0), Chern number: 2
# s1*s3 -> s1: (-1, 1, 0), Chern number: 2
# s1*s3 -> s1*s3*s2: (1, 0, 0), Chern number: 2
# s1*s3 -> s3: (0, -1, 1), Chern number: 0, con-as: ZZRingElem[0, -1, 2, -1]
# s2*s3 -> s2: (-1, 1, 0), Chern number: 2
# s2*s3 -> s2*s3*s2: (1, 0, 0), Chern number: 2
# s2*s3 -> s3: (0, 1, 0), Chern number: 4
# s1*s2*s3 -> s1*s2: (-1, 1, 0), Chern number: 2
# s1*s2*s3 -> s1*s2*s3*s2: (1, 0, 0), Chern number: 2
# s1*s2*s3 -> s1*s3: (0, 1, 0), Chern number: 4
# s1*s2*s3 -> s2*s3: (0, 0, 1), Chern number: 4

# (Root system of type A3, RootSpaceElem[], "s2*s1*s3"):
# Schubert variety is: GKM graph with 8 nodes and valency 3
# s1 -> id: (0, 0, 1), Chern number: 2
# s2*s1 -> s1: (-1, 1, 1), Chern number: 2
# s2 -> id: (-1, 1, 0), Chern number: 0, con-as: ZZRingElem[2, -1, -1]
# s2 -> s2*s1: (0, 0, 1), Chern number: 2
# s3 -> id: (1, 0, 0), Chern number: 2
# s1*s3 -> s1: (1, 0, 0), Chern number: 2
# s1*s3 -> s3: (0, 0, 1), Chern number: 2
# s2*s1*s3 -> s2*s1: (1, 0, 0), Chern number: 2
# s2*s1*s3 -> s1*s3: (0, 1, 1), Chern number: 4
# s2*s3 -> s2: (1, 0, 0), Chern number: 2
# s2*s3 -> s3: (0, 1, 0), Chern number: 2
# s2*s3 -> s2*s1*s3: (0, 0, 1), Chern number: 2

# (Root system of type A3, RootSpaceElem[a_1], "s1*s3*s2"):
# Schubert variety is: GKM graph with 6 nodes and valency 3
# s2 -> id: (0, 1), Chern number: 3
# s1*s2 -> id: (0, 1), Chern number: 3
# s1*s2 -> s2: (0, 1), Chern number: 3
# s3*s2 -> s2: (1, 1), Chern number: 3
# s1*s3*s2 -> s1*s2: (1, 1), Chern number: 3
# s1*s3*s2 -> s3*s2: (0, 1), Chern number: 3
# s3 -> id: (1, 0), Chern number: 0, con-as: ZZRingElem[2, -1, -1]
# s3 -> s3*s2: (0, 1), Chern number: 3
# s3 -> s1*s3*s2: (0, 1), Chern number: 3

# (Root system of type A3, RootSpaceElem[a_3], "s1*s3*s2"):
# Schubert variety is: GKM graph with 6 nodes and valency 3
# s1 -> id: (-1, 1), Chern number: 0, con-as: ZZRingElem[2, -1, -1]
# s2 -> id: (1, 0), Chern number: 3
# s1*s2 -> s1: (1, 0), Chern number: 3
# s1*s2 -> s2: (0, 1), Chern number: 3
# s3*s2 -> id: (1, 0), Chern number: 3
# s3*s2 -> s2: (1, 0), Chern number: 3
# s1*s3*s2 -> s1: (1, 0), Chern number: 3
# s1*s3*s2 -> s1*s2: (1, 0), Chern number: 3
# s1*s3*s2 -> s3*s2: (0, 1), Chern number: 3
# 5-element Vector{Any}:
#  (Root system of type A3, RootSpaceElem[], "s1*s3*s2*s1", 0)
#  (Root system of type A3, RootSpaceElem[], "s1*s2*s3*s2", 0)
#  (Root system of type A3, RootSpaceElem[], "s2*s1*s3", 0)
#  (Root system of type A3, RootSpaceElem[a_1], "s1*s3*s2", 0)
#  (Root system of type A3, RootSpaceElem[a_3], "s1*s3*s2", 0)

#  julia> search_nonpositive_Schubert(root_system(:G, 2); printFully=true)

# (Root system of type G2, RootSpaceElem[], "s1*s2"):
# Schubert variety is: GKM graph with 4 nodes and valency 2
# s1 -> id: (-3, 1), Chern number: -1
# s2 -> id: (1, 0), Chern number: 2
# s1*s2 -> s1: (1, 0), Chern number: 2
# s1*s2 -> s2: (0, 1), Chern number: 5

# (Root system of type G2, RootSpaceElem[], "s2*s1*s2"):
# Schubert variety is: GKM graph with 6 nodes and valency 3
# s1 -> id: (-2, 1), Chern number: 0, con-as: ZZRingElem[2, 0, -2]
# s2*s1 -> s1: (-1, 1), Chern number: 2
# s2 -> id: (1, 0), Chern number: 2
# s2 -> s2*s1: (-2, 1), Chern number: 0, con-as: ZZRingElem[0, 2, -2]
# s1*s2 -> s1: (1, 0), Chern number: 2
# s1*s2 -> s2: (1, 1), Chern number: 6
# s2*s1*s2 -> id: (1, 1), Chern number: 6
# s2*s1*s2 -> s2*s1: (1, 0), Chern number: 2
# s2*s1*s2 -> s1*s2: (0, 1), Chern number: 4
# 2-element Vector{Any}:
#  (Root system of type G2, RootSpaceElem[], "s1*s2", -1)
#  (Root system of type G2, RootSpaceElem[], "s2*s1*s2", 0)