using nauty_jll

# Determine C compiler, empty if not set
cc = Sys.which("gcc") !== nothing ? "gcc" : (Sys.which("clang") !== nothing ? "clang" : get(ENV, "CC", nothing))

# Determine platform and set compiler and flags
shared_flag = Sys.isunix() ? "-shared" : (Sys.isapple() ? "-dynamiclib" : "")

if cc === nothing || shared_flag == ""
  @warn "Skipping sizeaut build: no suitable C compiler (gcc or clang) found or unsupported OS ($(Sys.KERNEL) $(Sys.MACHINE))"
  return
end

@info "Using C compiler: $cc"
const sizeaut_path  = joinpath(@__DIR__, "sizeaut.c")
const libsizeaut    = joinpath(@__DIR__, "libsizeaut." * Base.Libc.Libdl.dlext)
const libnauty_path = nauty_jll.libnauty
const include_path  = joinpath(nauty_jll.artifact_dir, "include")

@info "Compiling wrapper library at $libsizeaut for sizeaut"
run(`$cc -fPIC $shared_flag -o $libsizeaut $sizeaut_path -I$include_path $libnauty_path -lm`)
@info "Compilation done"
