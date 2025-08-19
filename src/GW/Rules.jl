import Base: *, //, /, ^, +, -, inv, one, zero
# import Base: *

struct EquivariantClass
    rule::Expr
    func::Function
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

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function *(ec1::EquivariantClass, n::Number)::EquivariantClass

    rule = quote
        $(ec1.rule) * $(n)
    end

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function *(n::Number, ec1::EquivariantClass)::EquivariantClass

    return ec1 * n
end

### Division
function //(ec1::EquivariantClass, ec2::EquivariantClass)::EquivariantClass

    rule = quote
        $(ec1.rule) // $(ec2.rule)
    end

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function //(ec1::EquivariantClass, n::Number)::EquivariantClass

    rule = quote
        $(ec1.rule) // $(n)
    end

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

### Sums
function +(ec1::EquivariantClass, ec2::EquivariantClass)::EquivariantClass

    rule = quote
        $(ec1.rule) + $(ec2.rule)
    end

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

function +(ec1::EquivariantClass, n::Number)::EquivariantClass #it makes sene only if n==0

    rule = quote
        $(ec1.rule) + $(n)
    end

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end
function +(n::Number, ec1::EquivariantClass)::EquivariantClass #it makes sene only if n==0

    return ec1 + n
end
### Minus
function -(ec1::EquivariantClass, ec2::EquivariantClass)::EquivariantClass

    rule = quote
        $(ec1.rule) - $(ec2.rule)
    end

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end
function -(ec1::EquivariantClass, n::Number)::EquivariantClass #it makes sene only if n==0

    rule = quote
        $(ec1.rule) - $(n)
    end

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end
function -(n::Number, ec1::EquivariantClass)::EquivariantClass #it makes sene only if n==0

    rule = quote
        $(n) - $(ec1.rule)
    end

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end
### Exponent
function ^(ec1::EquivariantClass, n::Number)::EquivariantClass

    rule = quote
        $(ec1.rule)^$(n)
    end

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end

### Constants
function one(ec1::EquivariantClass)::EquivariantClass

    rule = quote
        1
    end

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end
function inv(ec1::EquivariantClass)::EquivariantClass

    rule = quote
        $(ec1.rule)^(-1)
    end

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end
function zero(ec1::EquivariantClass)::EquivariantClass

    rule = quote
        0
    end

    return EquivariantClass(rule, eval(:((dt) -> $rule)))
end