##############################
### Equivariant classes ###
##############################

struct EquivariantClass
  rule::Expr
  func::Function
  # has_psi::Int64 # has_psi: 0-no, 1-yes, 2-sum_found, 3-prod_found
  sum_psi::Bool # true if a class with psi and a class without are summed
  prod_psi::Bool # true if two classes with psi are multiplied
  psi_exp::Vector{Int64} # exponents of psi classes at markings
  lambda_coef::Vector{Int64} # coefficients of lambda classes
end

# by default no psi classes and no lambda classes
function EquivariantClass(rule, func)::EquivariantClass
  EquivariantClass(rule, func, false, false, Int64[], Int64[])
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

  #### has_psi is used to check the validity of the operations
  
  sum_psi = ec1.sum_psi || ec2.sum_psi
  prod_psi = ec1.prod_psi || ec2.prod_psi || ( (!isempty(ec1.psi_exp)) && (!isempty(ec2.psi_exp)) )
  psi_exp = isempty(ec1.psi_exp) ? ec2.psi_exp : ec1.psi_exp
  lambda_coef = isempty(ec1.lambda_coef) ? ec2.lambda_coef : ec1.lambda_coef

  return EquivariantClass(rule, eval(:((dt) -> $rule)), sum_psi, prod_psi, psi_exp, lambda_coef)
  # return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function *(ec1::EquivariantClass, n::Number)::EquivariantClass
  rule = quote
    $(ec1.rule) * $(n)
  end

  return EquivariantClass(rule, eval(:((dt) -> $rule)), ec1.sum_psi, ec1.prod_psi, ec1.psi_exp, ec1.lambda_coef)
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

  return EquivariantClass(rule, eval(:((dt) -> $rule)), ec1.sum_psi, ec1.prod_psi, ec1.psi_exp, ec1.lambda_coef)
end

### Sums
function +(ec1::EquivariantClass, ec2::EquivariantClass)::EquivariantClass
  rule = quote
    $(ec1.rule) + $(ec2.rule)
  end

  sum_psi = ec1.sum_psi || ec2.sum_psi || ( (isempty(ec1.psi_exp) != isempty(ec2.psi_exp)) )
  prod_psi = ec1.prod_psi || ec2.prod_psi

  return EquivariantClass(rule, eval(:((dt) -> $rule)), sum_psi, prod_psi, ec1.psi_exp, ec1.lambda_coef)
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
  # rule = quote
  #   $(ec1.rule) - $(ec2.rule)
  # end

  # has_psi = ec1.has_psi == ec2.has_psi ? ec1.has_psi : max(ec1.has_psi, ec2.has_psi, 2)

  # return EquivariantClass(rule, eval(:((dt) -> $rule)), has_psi)
  return ec1 + ( -1 * ec2 )
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

  prod_psi = ec1.prod_psi || ( (!isempty(ec1.psi_exp)) && (n != 0) )
  return EquivariantClass(rule, eval(:((dt) -> $rule)), ec1.sum_psi, prod_psi, ec1.psi_exp, ec1.lambda_coef)
end

### Constants
function one(ec1::EquivariantClass)::EquivariantClass
  rule = quote
    1
  end

  return EquivariantClass(rule, eval(:((dt) -> $rule)), false, false, Int64[], Int64[])
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

  return EquivariantClass(rule, eval(:((dt) -> $rule)), false, false, Int64[], Int64[])
end