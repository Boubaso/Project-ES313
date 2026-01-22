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
println("Starting  simulation...")
println("------------------------------------------------")
println("------------------------------------------------")

# 2D Frank–Kamenetskii Reaction–Diffusion Model (Julia)
# Graphical visualization using Plots.jl

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

# Numerical grid
    Nx, Ny = 201, 201
    Δx = 1e-3            # m
    Δt = 1e-1            # s
    tmax = 15             # s
    estimated_iterations = Int(tmax/Δt)
    println("Estimated max iterations: $estimated_iterations")

    geometry = :plane  # :plane or :cylinder
    simulation_video_time = 5.0  # seconds of simulation video
    fps = 10 

    β = material.α * Δt / Δx^2
    println("Diffusion number β = $β (stable if ≤ 0.05)")
    if β > 0.05
        @error "Unstable simulation parameters: reduce Δt or increase Δx"
    end

# Initial condition: localized hot spot
    T_ambient = 300.0    # K
    T_hot = 530.0        # K  (adjust to see sub/supercritical behavior)
    R_hot = 1e-3         # m, radius circle or radius of line (diameter/2)

    T = fill(T_ambient, Nx, Ny)
    Tnew = similar(T)
    cx, cy = div(Nx,2), div(Ny,2)

    initial_sit = "line"  # "line" or "circle"
        # Initial situation: "circle" for circular hot spot
            # considered an impuls in the center, ESD 
        # Initial situation: "line" for line hot spot
            # considered a wire that gives a voltage input as heat



#############################Separation input and code########################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
    




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

# Plot setup
plt_heat = heatmap(
    (1:Ny) .* Δx,
    (1:Nx) .* Δx,
    T',
    aspect_ratio=1,
    c=:inferno,
    clims=(T_ambient,T_hot*2),
    title="Temperature field (K)",
    xlabel="x ($(Δx)m)",
    ylabel="y ($(Δx)m)"
)

plt_graph = plot(
    xlabel="Time (s)",
    ylabel="Max temperature (K)",
    title="Ignition monitoring",
    legend=false
)

full_plot = plot(plt_heat, plt_graph, layout=(1,2), size=(1000,650))
display(full_plot)

###### FUNCTIONS ######
@inline function reaction(T, material)
    return (material.Q*material.A/(material.ρ*material.cp)) * exp(-material.E/(material.Rg*T))
end

@inline function laplacian(T,i,j;Δx=Δx)
    return (T[i+1,j] + T[i-1,j] + T[i,j+1] + T[i,j-1] - 4T[i,j]) / Δx^2 # FDM 5-point stencil
end

function output_record(fire, iter, t, Tmax, Tavg, Δx, Δt, β, initial_sit, geometry, material, sim_time)
    # Create detailed simulation record
            record = Dict(
                "result" => fire ? "FIRE" : "NO FIRE",
                "final_iteration" => iter,
                "final_time" => t,
                "max_temperature" => Tmax,
                "avg_temperature" => Tavg,
                "grid_size" => (Nx, Ny),
                "spatial_step" => Δx,
                "time_step" => Δt,
                "diffusion_number" => β,
                "initial_condition" => initial_sit,
                "geometry" => geometry,
                "material" => material,
                "simulation_time" => sim_time
            )
            println("\n========== SIMULATION SUMMARY ==========")
            println("Outcome: $(record["result"])")
            println("Final iteration: $(record["final_iteration"])")
            println("Final time: $(round(record["final_time"], digits=3)) s")
            println("Max temperature end: $(round(record["max_temperature"], digits=2)) K")
            println("Avg temperature end: $(round(record["avg_temperature"], digits=2)) K")
            println("Grid size: $(record["grid_size"])")
            println("Spatial step: $(record["spatial_step"]) m")
            println("Time step: $(record["time_step"]) s")
            println("Diffusion number β: $(round(record["diffusion_number"], digits=4))")
            println("Simulation time: $(record["simulation_time"]) s")
            println("========================================\n")

            # Save record to text file
            open("ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT/CA numerical_output.txt", "w") do io
                println(io, "========== SIMULATION SUMMARY ==========")
                println(io, "Outcome: $(record["result"])")
                println(io, "Final iteration: $(record["final_iteration"])")
                println(io, "Final time: $(round(record["final_time"], digits=3)) s")
                println(io, "Max temperature end: $(round(record["max_temperature"], digits=2)) K")
                println(io, "Avg temperature end: $(round(record["avg_temperature"], digits=2)) K")
                println(io, "Grid size: $(record["grid_size"])")
                println(io, "Spatial step: $(record["spatial_step"]) m")
                println(io, "Time step: $(record["time_step"]) s")
                println(io, "Diffusion number β: $(round(record["diffusion_number"], digits=4))")
                println(io, "Initial condition: $(record["initial_condition"])")
                println(io, "Geometry: $(record["geometry"])")
                println(io, "Simulation time: $(record["simulation_time"]) s")
                println(io, "========================================")
            end

end


###### TIME INTEGRATION ######

global t = 0.0

iter = Int(1)
times = Float64[0.0]
Tmax_hist = Float64[maximum(T)]
Tavg_hist = Float64[mean(T)]
Tmax_dt_hist = Float64[0.0]

FULLPLOT = [full_plot]

println("Starting time integration...")
t_start = Base.time()

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

    global T .= Tnew
    global t += Δt
    global iter += 1

    local Tmax = maximum(T)
    local Tavg = mean(T)
    local Tmax_dt = (Tmax - Tmax_hist[end]) / Δt
    
    push!(times, t)
    println("iteration $(iter)/ $(estimated_iterations): t=$(round(t,digits=3))s")
    push!(Tmax_hist, Tmax)
    push!(Tavg_hist, Tavg)
    push!(Tmax_dt_hist, Tmax_dt)

    if iter % round(estimated_iterations / (simulation_video_time * fps)) == 0
        global plt_heat = heatmap!(plt_heat, T', clims=(300,1200), title="Temperature field (K) - t=$(round(t,digits=2))s")
        global plt_graph = plot!(plt_graph, times, Tmax_hist)
        global full_plot = plot(plt_heat, plt_graph, layout=(1,2), size=(1000,650))
        push!(FULLPLOT, full_plot)
        println("t = $(round(t,digits=3)) s, Tmax = $(round(Tmax,digits=2)) K, ΔTmax = $(round(Tmax_dt,digits=4)) K")
        display(full_plot)
    end

    if iter == 178
        @info "Vive la 178 !!"
    end

    # Stop criteria
    if t > tmax
        local t_end = Base.time()
        @info "❄️ NO FIRE and tmax reached at t = $(round(t,digits=3)) s, Tmax = $(round(Tmax,digits=2)) K"
        output_record(false,iter, t, Tmax, Tavg, Δx, Δt, β, initial_sit, geometry, material, (round(t_end - t_start, digits=2)))
        break
    elseif (Tmax_dt - Tmax_dt_hist[end-1]) < -20 && Tmax < 3000
        local t_end = time()
        @info "❄️ NO FIRE + stability at t = $(round(t,digits=3)) s, acc Tmax < 20 K/s"
        output_record(false, iter, t, Tmax, Tavg, Δx, Δt, β, initial_sit, geometry, material, (round(t_end - t_start, digits=2)))
        break
    elseif (Tmax_dt - Tmax_dt_hist[end-1]) > 20 # K/s
        t_end = time()
        @info "🔥 FIRE detected at t = $(round(t,digits=3)) s, acc Tmax > 20 K/s"
        output_record(true, iter, t, Tmax, Tavg, Δx, Δt, β, initial_sit, geometry, material, (round(t_end - t_start, digits=2)))
        break
    elseif Tmax > 10^4
        local t_end = time()
        @info "🔥 FIRE detected at t = $(round(t,digits=3)) s, Tmax too big"
        output_record(true, iter, t, Tmax, Tavg, Δx, Δt, β, initial_sit, geometry, material, (round(t_end - t_start, digits=2)))
        break
    elseif iter > estimated_iterations * 1.1
        local t_end = time()
        @warn "Stopping simulation: exceeded maximum expected iterations"
        output_record(false, iter, t, Tmax, Tavg, Δx, Δt, β, initial_sit, geometry, material, (round(t_end - t_start, digits=2)))
        break
    end
end


#############################Separation output and code#######################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

# Save final plots
savefig(FULLPLOT[end], "ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT/CA numerical output.png")
# Save animation
function save_animation(Figs, Dirtry = "animation.mp4", Dircatch = "animation.gif")
    anim = Plots.Animation()
    for p in Figs
        Plots.frame(anim, p)
    end
    try
        Plots.mp4(anim, Dirtry, fps=10)
    catch
        Plots.gif(anim, Dircatch, fps=10)
    end
end

if !isempty(FULLPLOT)      
save_animation(FULLPLOT, "ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT/CA numerical output.mp4", "ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT/CA numerical output.gif")
end  

println("Code completed.")


