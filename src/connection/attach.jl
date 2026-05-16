function set_connection!(g::GKMGraph, new::AbstractGKMConnection{R}) where {R}
  old = g.connection

  # If we already have a canonical one, do nothing
  if old isa CanonicalConnection{R}
    return nothing
  end

  # If new is canonical OR nothing stored yet, accept it
  if new isa CanonicalConnection{R} || isnothing(old)
    g.connection = new
  end

  return nothing
end

function get_connection(g::GKMGraph; require_canonical::Bool=false)
  if isnothing(g.connection)
    # compute fallback first
    g.connection = compute_algorithmic_connection(g)
  end

  if require_canonical && !(g.connection isa CanonicalConnection)
    c = try_compute_canonical_connection(g)
    if c !== nothing
      g.connection = c
    end
  end

  return g.connection
end