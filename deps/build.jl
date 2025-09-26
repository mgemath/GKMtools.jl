using nauty_jll

const sizeaut_path = joinpath(@__DIR__, "sizeaut.c")
const libsizeaut = joinpath(@__DIR__, "libsizeaut." * Base.Libc.Libdl.dlext)
const libnauty_path = nauty_jll.libnauty
const include_path = joinpath(nauty_jll.artifact_dir, "include")

@info "Compiling wrapper library at $libsizeaut for sizeaut"
run(`gcc -fPIC -shared -o $libsizeaut $sizeaut_path -I$include_path $libnauty_path -lm`)

@info "Compilation done"
