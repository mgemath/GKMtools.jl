##############################
### Equivariant classes ###
##############################

struct EquivariantClass
  rule::Expr
  func::Function
  has_psi::Int64 # has_psi: 0-no, 1-yes, 2-sum_found, 3-prod_found 
end

# by default has_psi = 0, the only exception is the psi class
function EquivariantClass(rule, func)::EquivariantClass
  EquivariantClass(rule, func, 0)
end

function Base.show(::IO, ::EquivariantClass) end

#########################################
### Operations of equivariant classes ###
#########################################

### Products
function *(ec1::EquivariantClass, ec2::EquivariantClass)::EquivariantClass
  rule = quote
    $(ec1.rule) * $(ec2.rule)
  end

  has_psi = (ec1.has_psi==1 && ec2.has_psi==1) ? 3 : max(ec1.has_psi, ec2.has_psi)

  return EquivariantClass(rule, eval(:((dt) -> $rule)), has_psi)
  # return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function *(ec1::EquivariantClass, n::Number)::EquivariantClass
  rule = quote
    $(ec1.rule) * $(n)
  end

  return EquivariantClass(rule, eval(:((dt) -> $rule)), ec1.has_psi)
  # return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function *(n::Number, ec1::EquivariantClass)::EquivariantClass
  return ec1 * n
end

### Division
# function //(ec1::EquivariantClass, ec2::EquivariantClass)::EquivariantClass
#   rule = quote
#     $(ec1.rule)//$(ec2.rule)
#   end

#   return EquivariantClass(rule, eval(:((dt) -> $rule)))
# end

function //(ec1::EquivariantClass, n::Number)::EquivariantClass
  rule = quote
    $(ec1.rule)//$(n)
  end

  return EquivariantClass(rule, eval(:((dt) -> $rule)), ec1.has_psi)
end

### Sums
function +(ec1::EquivariantClass, ec2::EquivariantClass)::EquivariantClass
  rule = quote
    $(ec1.rule) + $(ec2.rule)
  end

  has_psi = ec1.has_psi == ec2.has_psi ? ec1.has_psi : max(ec1.has_psi, ec2.has_psi, 2)

  return EquivariantClass(rule, eval(:((dt) -> $rule)), has_psi)
end

# function +(ec1::EquivariantClass, n::Number)::EquivariantClass #it makes sene only if n==0
#   rule = quote
#     $(ec1.rule) + $(n)
#   end

#   return EquivariantClass(rule, eval(:((dt) -> $rule)))
# end
# function +(n::Number, ec1::EquivariantClass)::EquivariantClass #it makes sene only if n==0
#   return ec1 + n
# end
### Minus
function -(ec1::EquivariantClass, ec2::EquivariantClass)::EquivariantClass
  rule = quote
    $(ec1.rule) - $(ec2.rule)
  end

  has_psi = ec1.has_psi == ec2.has_psi ? ec1.has_psi : max(ec1.has_psi, ec2.has_psi, 2)

  return EquivariantClass(rule, eval(:((dt) -> $rule)), has_psi)
end
# function -(ec1::EquivariantClass, n::Number)::EquivariantClass #it makes sene only if n==0
#   rule = quote
#     $(ec1.rule) - $(n)
#   end

#   return EquivariantClass(rule, eval(:((dt) -> $rule)))
# end
# function -(n::Number, ec1::EquivariantClass)::EquivariantClass #it makes sene only if n==0
#   rule = quote
#     $(n) - $(ec1.rule)
#   end

#   return EquivariantClass(rule, eval(:((dt) -> $rule)))
# end
### Exponent
function ^(ec1::EquivariantClass, n::Number)::EquivariantClass
  rule = quote
    $(ec1.rule)^$(n)
  end

  return EquivariantClass(rule, eval(:((dt) -> $rule)), ec1.has_psi)
end

### Constants
function one(ec1::EquivariantClass)::EquivariantClass
  rule = quote
    1
  end

  return EquivariantClass(rule, eval(:((dt) -> $rule)), 0)
end
# function inv(ec1::EquivariantClass)::EquivariantClass
#   rule = quote
#     $(ec1.rule)^(-1)
#   end

#   return EquivariantClass(rule, eval(:((dt) -> $rule)))
# end
function zero(ec1::EquivariantClass)::EquivariantClass
  rule = quote
    0
  end

  return EquivariantClass(rule, eval(:((dt) -> $rule)), 0)
end