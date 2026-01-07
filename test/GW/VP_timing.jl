n_repeats = 1

F = gkm_3d_twisted_flag();
beta = curve_class(F, "3", "4");

@time begin
    for _ in 1:n_repeats
        for d in 1:3
            for g in 1:2
                gw = gromov_witten(F, d*beta, 0, class_one(); g=g, show_bar=false)
                println("Genus $g, degree $d: $gw")
            end
        end
    end
end