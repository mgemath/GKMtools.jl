# This file is part of GKMtools.jl, licensed under the MIT License (MIT).

## GKM graphs and related types
abstract type AbstractGKMGraph{R, V, F} end
abstract type AbstractStackyGKMGraph{R, V, F} <: AbstractGKMGraph{R, V, F} end
abstract type AbstractOrbifoldGKMGraph{R, V, F} <: AbstractStackyGKMGraph{R, V, F} end

# Vertex and flag types for GKM graphs
abstract type AbstractVertex end

abstract type AbstractFlagWeight{R} end
abstract type AbstractOrbifoldFlagWeight{R} <: AbstractFlagWeight{R} end

abstract type AbstractGKMConnection{R} end

# Isotropy data for vertices and edge multiplicities for stacky GKM graphs
abstract type AbstractIsotropy{Gr, Rep} end
abstract type AbstractVertexIsotropy{Gr_type, Rep} <: AbstractIsotropy{Gr_type, Rep} end
abstract type AbstractFlagIsotropy{Gr_type, Rep, Embedding} <: AbstractIsotropy{Gr_type, Rep} end

# Examples of stacky fans and their associated GKM graphs
abstract type AbstractStackyFan end
abstract type AbstractStackyCone end