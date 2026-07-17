function _validate_homogeneous_connection_option(connection::Symbol)
  connection in (:geometric, :combinatorial) || throw(
    ArgumentError(
      "connection must be :geometric or :combinatorial; got $(repr(connection))",
    ),
  )
  return connection
end
