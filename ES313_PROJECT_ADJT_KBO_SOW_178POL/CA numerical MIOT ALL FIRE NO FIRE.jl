try 
    using Plots
    using Statistics
catch
    import Pkg; Pkg.add("Plots")
    import Pkg; Pkg.add("Statistics")
    using Plots
    using Statistics
end
gr()

println("Packages loaded successfully.")
println("------------------------------------------------")
println("------------------------------------------------")

###### INPUT ######
# material chemistry
    struct Material
        ρ::Float64      # kg/m^3
        cp::Float64       # J/(kg·K)
        k::Float64        # W/(m·K)
        α::Float64    # thermal diffusivity (m^2/s)
        E::Float64        # J/mol
        A::Float64        # 1/s
        Q::Float64        # J/kg
        Rg::Float64       # J/(mol·K)
    end

    const RDX = Material(
        1800.0,           # ρ
        1250.0,           # cp
        0.5,              # k
        0.5 / (1800.0 * 1250.0),  # α
        190e3,            # E
        1e18,             # A
        5e6,              # Q
        8.314             # Rg
    )

    material = RDX
    material_name = "RDX"

# Numerical grid
    Nx, Ny = 201, 201
    Δx = 1e-3            # m
    Δt = 1e-2            # s
    tmax = 30             # s
    estimated_iterations = Int(tmax/Δt)
    geometry = :plane  # :plane or :cylinder
    β = material.α * Δt / Δx^2
    if β > 0.05
        @error "Unstable simulation parameters: reduce Δt or increase Δx"
    end

# Initial condition: localized hot spot
    T_ambient = 300.0    # K
    R_hot = 5e-4         # m, radius circle or radius of line (diameter/2)
    initial_sit = "line"  # "line" or "circle"

# Data storage for MIOT
    IDX = []
    TEMPS = []
    FIRE = []
    TIMES = []

#############################Separation input and code########################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

###### FUNCTIONS ######
@inline function reaction(T, material)
    return (material.Q*material.A/(material.ρ*material.cp)) * exp(-material.E/(material.Rg*T))
end

@inline function laplacian(T,i,j;Δx=Δx)
    return (T[i+1,j] + T[i-1,j] + T[i,j+1] + T[i,j-1] - 4T[i,j]) / Δx^2 # FDM 5-point stencil
end

@inline function initiation(Nx, Ny, R_hot, T_hot, initial_sit; Δx=Δx)
    T = fill(T_ambient, Nx, Ny)
    cx, cy = div(Nx,2), div(Ny,2)
    if initial_sit == "circle"
            for i in 1:Nx, j in 1:Ny
            local r = sqrt(((i-cx)*Δx)^2 + ((j-cy)*Δx)^2)
                if r ≤ R_hot
                    T[i,j] = T_hot
                end
            end
        elseif initial_sit == "line"
            for i in 1:Nx, j in 1:Ny
            local x = (i - cx)*Δx
                if abs(x) <= R_hot
                    T[i,j] = T_hot
                end
            end
        else
            error("Invalid initial situation, must be 'circle' or 'line'")
    end
    return T
end

@inline function simulation(T, Nx, Ny, Δx, Δt, tmax, estimated_iterations, geometry, material, R_hot, T_hot, initial_sit)
    t = 0.0
    iter = Int(1)
    times = Float64[0.0]
    Tmax_hist = Float64[maximum(T)]
    Tmax_dt_hist = Float64[0.0]
   
    Tnew = similar(T)
    Fire = false

    while t < 1.1*tmax
        for i in 2:Nx-1, j in 2:Ny-1
            Tnew[i,j] =
                T[i,j] +
                Δt * (material.α * laplacian(T,i,j,Δx=Δx) +
                reaction(T[i,j], material)
                )
        end
        # to be considered infinite in approximation    
            Tnew[:,1]   .= Tnew[:,2]
            Tnew[:,end] .= Tnew[:,end-1]
            Tnew[1,:]   .= Tnew[2,:] 
        if geometry == :plane
            Tnew[end,:] .= Tnew[end-1,:]
        elseif geometry == :cylinder
            Tnew[end,:] .= Tnew[1,:]  # wrap-around in θ
        else
            error("Invalid geometry, must be :plane or :cylinder")
        end
        
        T .= Tnew
        t += Δt
        iter += 1
        Tmax = maximum(T)
        Tmax_dt = (Tmax - Tmax_hist[end]) / Δt
        
        push!(times, t)
        push!(Tmax_hist, Tmax)
        push!(Tmax_dt_hist, Tmax_dt)

        # Stop criteria
        if t > tmax
            break
        elseif (Tmax_dt - Tmax_dt_hist[end-1]) < -5 && Tmax < 1500
            break
        elseif (Tmax_dt - Tmax_dt_hist[end-1]) > 20 || Tmax > 1500
            Fire = true
            break
        elseif iter > estimated_iterations * 1.1
            break
        end
    end
    return Fire
end

function fires(T_hot, Nx, Ny, Δx, Δt, tmax, estimated_iterations, geometry, material, R_hot, initial_sit)
    T = initiation(Nx, Ny, R_hot, T_hot, initial_sit)
    return simulation(
        T, Nx, Ny, Δx, Δt, tmax, estimated_iterations,
        geometry, material, R_hot, T_hot, initial_sit
    )
end

function find_critical_temperature_with_tracking(Tmin, Tmax; Tol=1e-2, maxiter=50)
        fmin = fires(Tmin, Nx, Ny, Δx, Δt, tmax, estimated_iterations, geometry, material, R_hot, initial_sit)
        fmax = fires(Tmax, Nx, Ny, Δx, Δt, tmax, estimated_iterations, geometry, material, R_hot, initial_sit)

        if fmin
            error("Already fires at Tmin = $Tmin")
        end
        if !fmax
            error("Does not fire at Tmax = $Tmax")
        end

        iter = 0
        while (Tmax - Tmin) > Tol && iter < maxiter
            Tmid = 0.5 * (Tmin + Tmax)
            
            # Time the simulation
            start_time = Base.time()
            result = fires(Tmid, Nx, Ny, Δx, Δt, tmax, estimated_iterations, geometry, material, R_hot, initial_sit)
            elapsed_time = Base.time() - start_time
            
            # Store results
            push!(IDX, iter + 1)
            push!(TEMPS, Tmid)
            push!(FIRE, result)
            push!(TIMES, elapsed_time)

            if result
                Tmax = Tmid   # ignition → lower it
            else
                Tmin = Tmid   # no ignition → raise it
            end

            iter += 1
        end

        return 0.5 * (Tmin + Tmax)
end

function format_number(x)
    replace(string(x), "." => ",")
end
########### MIOT (multiple iterations, optimal temperature #####################

Tol = 1.0e-2
T_hot_min = 350.0    # K  (adjust to see sub/supercritical behavior)
T_hot_max = 800.0    # K

Tcrit = find_critical_temperature_with_tracking(T_hot_min, T_hot_max; Tol=Tol)
try 
    using DataFrames
    using CSV
catch
    import Pkg; Pkg.add("DataFrames")
    import Pkg; Pkg.add("CSV")
end

notes = ["Nx=$Nx, Ny=$Ny, dx=$Δx m, dt=$Δt s, tmax=$tmax s, geometry=$geometry, initial_sit=$initial_sit, R_hot=$R_hot m, material =$material_name"]
NOTES = [i == 1 ? notes[1] : "" for i in 1:length(IDX)]
df = DataFrame(
        idx = [format_number(p) for p in IDX],
        T_hot_K = [format_number(p) for p in TEMPS],
        Fire = [format_number(p) for p in FIRE], 
        Time_s = [format_number(p) for p in TIMES],
        Note = NOTES
    )
CSV.write("ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT/CA MIOT output.csv", df; delim=';')

println("------------------------------------------------")
println("Critical ignition temperature:")
println("T_hot,crit ≈ $(round(Tcrit, digits=3)) K")
println("after $(length(IDX)) iterations")
println("------------------------------------------------")
 

println("Code completed.")