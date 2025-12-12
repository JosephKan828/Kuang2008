module IO

using HDF5
export load_background, load_vertical_modes, load_domain, load_wavenumbers, save_state

function load_background(path::String)
    h5open(path, "r") do f
        ρ0 = read(f, "ρ0")
        p0 = read(f, "p0")
        T0 = read(f, "T0")
        z  = read(f, "z")
        return ρ0, p0, T0, z
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

function load_wavenumbers(path::String)
    h5open(path, "r") do f
        k = read(f, "wavenumber")
        return k
    end
end

function save_state(path::String, state, t, k, varnames)
    h5open(path, "w") do f
        write(f, "state", state)
        write(f, "time", t)
        write(f, "wavenumber", k)
        write(f, "variables", varnames)

        attributes(f["time"])["units"] = "day"
        attributes(f["wavenumber"])["units"] = "None"
    end
end

end # module
