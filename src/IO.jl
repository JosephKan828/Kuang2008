module IO

using HDF5
export load_background, load_vertical_modes, load_domain, load_inv_mat, save_state, save_optrs

function load_background(path::String)
    h5open(path, "r") do f
        ρ0 = read(f, "rho_bar")
        p0 = read(f, "pbar")
        T0 = read(f, "Tbar")
        return ρ0, p0, T0
    end
end

function load_vertical_modes(path::String)
    h5open(path, "r") do f
        G1 = read(f, "G1")
        G2 = read(f, "G2")
        return G1, G2
    end
end

function load_domain(path::String)
    h5open(path, "r") do f
        x = read(f, "x")
        z = read(f, "z")
        t = read(f, "t")
        return x, z, t
    end
end

function load_inv_mat(path::String)
    h5open(path, "r") do f
        λ = read(f, "lambda")
        kcal = read(f, "wnum_cal")
        kdis = read(f, "wnum_dis")
        Finv = read(f, "F_inv")
        return λ, kcal, kdis, Finv
    end
end

function save_state(path::String, state, t, k, varnames)
    println("Saving state to ", path)
    h5open(path, "w") do f
        write(f, "state", state)
        write(f, "time", t)
        write(f, "wavenumber", k)
        write(f, "variables", varnames)

        attributes(f["time"])["units"] = "day"
        attributes(f["wavenumber"])["units"] = "None"
    end
end

function save_optrs( path::String, optrs, k )
    println( "Saving operators to ", path )

    h5open( path, "w" ) do f
        write( f, "operators", optrs )
        write( f, "wavenumber", k )
    end
end

end # module
