try
    using Plots
    using Statistics
catch
    import Pkg; Pkg.add("Plots"); using Plots
    import Pkg; Pkg.add("Statistics"); using Statistics
end
# INPUT
    #material properties (COPPER)
    α = 111*10^(-6)     # thermal diffusivity m^2/s
    L = 1.0             # length of the rod [m]
    nx = 65             # number of knodes
    t_final = 600
    speed = 10       # speed factor for simulation
# cODE
fps = 10
dx = L/(nx-1)
dt = 1/fps
β = α * dt / (dx^2)
if β >= 0.05 error("β too large for stability condition; increase nx or decrease fps") end
nt = t_final*fps
U = zeros(Float64, nt, nx)
U[1, :] .= 298.15
U[:, 1] .= 1000
U[:, end] .= 1000
# Preallocate array for average temperatures
avg_temps = Float64[]
time_points = Float64[]

for t in 1:(nt-1)
    for i in 2:(nx-1)
        U[t+1, i] = U[t, i] + β * (U[t, i+1] - 2*U[t, i] + U[t, i-1])
    end
    # Calculate average temperature at speed intervals
    if t % speed == 0
        push!(avg_temps, mean(U[t, :]))
        push!(time_points, (t-1)*dt)
    end
end

x = range(0, L, length=nx)
clims = (minimum(U), maximum(U))
anim = @animate for t in 1:speed:nt
    z = repeat(reshape(U[t, :], 1, nx), 2, 1)
    heatmap(x, [0,1], z; 
        xlabel = "x [m]", yticks = false, xlim = (0, L), title = "t = $(round((t-1)*dt, digits=2)) s",
        aspect_ratio = 0.15, size = (800, 200), 
        colorbar = true, colorbar_title = "Temperature [K]", clims = clims, c=:thermal
    )
end

# Plot average temperature over time
plot(time_points, avg_temps, 
    xlabel="Time [s]", ylabel="Average Temperature [K]",
    title="Average Temperature vs Time",
    linewidth=2, legend=false,
    size=(600,400))
savefig("ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT\\SOLN_1D_Conduction.png")
mp4(anim, "ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT\\SOLN_1D_Conduction.mp4", fps = fps)