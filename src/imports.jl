import Base: *, //, /, ^, +, -, inv, one, zero, length, <, >, isless
import Oscar.direct_sum, Oscar.line_bundle, Oscar.tensor_product
import Oscar.IntersectionTheory.tangent_bundle
import Oscar.IntersectionTheory.cotangent_bundle
import Oscar.total_space
import Oscar: projective_space
import Oscar: isvalid, blow_up, point_class, integrate, is_effective, chern_number, betti_numbers, rank, dual, projectivization, tangent_bundle, cotangent_bundle, chern_class, chern_classes, det
import Oscar: schubert_class, schubert_classes

# Graph-functions
import Oscar: Graph, Edge, all_neighbors, src, dst, add_vertex!, add_edge!, isvalid, is_connected, is_simple, is_loopless, neighbors, degree, indegree, outdegree, has_edge, has_vertex, vertices, edges, nv, ne

# Cone and toric variety functions
import Oscar: cones, maximal_cones, rays, dim, n_rays, polarize