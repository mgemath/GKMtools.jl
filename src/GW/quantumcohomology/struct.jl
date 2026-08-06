@doc raw"""
    SmallEquivariantQuantumCohomology

Parent type for a finite Novikov truncation of small equivariant quantum
cohomology. Use [`small_equivariant_quantum_cohomology`](@ref) to construct it.
"""
struct SmallEquivariantQuantumCohomology{G}
  graph::G
  degrees::Vector{CurveClass}
end

@doc raw"""
    SmallEquivariantQuantumClass

Element type for a finite sum ``\sum_\beta a_\beta q^\beta`` with localized
[`GKMClass`](@ref) coefficients. Use [`quantum_class`](@ref) to construct it.
"""
struct SmallEquivariantQuantumClass{P}
  parent::P
  coefficients::Dict{CurveClass,GKMClass}
end

Base.parent(c::SmallEquivariantQuantumClass) = c.parent
graph(QH::SmallEquivariantQuantumCohomology) = QH.graph
graph(c::SmallEquivariantQuantumClass) = graph(parent(c))

@doc raw"""
    truncation_degrees(QH)
    truncation_degrees(c)

Return a copy of the Novikov degrees retained by a quantum-cohomology parent
or quantum class.
"""
truncation_degrees(QH::SmallEquivariantQuantumCohomology) = copy(QH.degrees)
truncation_degrees(c::SmallEquivariantQuantumClass) = truncation_degrees(parent(c))

function Base.show(io::IO, QH::SmallEquivariantQuantumCohomology)
  print(io, "small equivariant quantum cohomology truncated to ",
    length(QH.degrees), " Novikov degrees")
end

function Base.show(io::IO, c::SmallEquivariantQuantumClass)
  isempty(c.coefficients) && return print(io, "0")
  terms = sort!(collect(c.coefficients); by=x -> string(first(x)))
  for (i, (beta, value)) in enumerate(terms)
    i > 1 && print(io, Oscar.is_terse(io) ? " + " : "\n + ")
    print(io, "(", value, ") q^(", beta, ")")
  end
end
Base.show(io::IO, ::MIME"text/plain", c::SmallEquivariantQuantumClass) = show(io, c)
