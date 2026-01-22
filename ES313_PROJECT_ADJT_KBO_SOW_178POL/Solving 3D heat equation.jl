try
    using Plots
    using LinearAlgebra
    using Statistics
catch
    import Pkg; Pkg.add("Plots")
    import Pkg; Pkg.add("LinearAlgebra")
    import Pkg; Pkg.add("Statistics")
    using Plots
    using LinearAlgebra
    using Statistics
end
# INPUT

α = 111*10^(-6)     # thermal diffusivity m^2/s
h = 10.0            # convective heat transfer coefficient W/m^2K
T_air = 298.15     # ambient temperature K

L = 1.0             # length of the surf [m]
nx = 20           # number of knodes
t_final = 600
speed = 100      # speed factor for simulation

fps = 10
dx = L/(nx-1)
dt = 1/fps
β = α * dt / (dx^2)
if β >= 0.05 error("β too large for stability condition; increase nx or decrease fps") end
nt = t_final*fps
U_vec = []
U = ones(Float64, nt, nx, nx, nx)*298.15
U[:, 1, :, :] .= 1000
U[:, end, :, :] .= 1000

#CONDUCTION
function conduction_euler_progressive(U, β, nx)
    U_new = copy(U)
    for i in 2:(nx-1), j in 2:(nx-1), k in 2:(nx-1)
        U_new[i, j, k] = U[i, j, k] + β * (U[i+1, j, k] + U[i-1, j, k] + U[i, j+1, k] + U[i, j-1, k] + U[i, j, k+1] + U[i, j, k-1] - 6*U[i, j, k])
    end
    return U_new
end

#CONVECTION
function convection_euler_progressive(U, h, dt, T_air, nx)
    U_new = copy(U)
    for i in nx, j in (1,nx), k in (1,nx)
        U_new[i, j, k] = U[i, j, k] + h * dt *(T_air - U[i, j, k])
    end
    return U_new
end

#HEAT TRANSFER (CONDUCTION + CONVECTION)
function heat_transfer(U, h, β, nx)
     U_new = copy(U)
     U_new = conduction_euler_progressive(U, β, nx)
     U_new = convection_euler_progressive(U_new, h, dt, T_air, nx)
    return U_new
end


#SIMULATION
avg_temp = []
for t in 1:(nt-1)
    U[t+1, :, :, :] = heat_transfer(U[t, :, :, :], h, β, nx)
    push!(avg_temp, mean(U[t+1, :, :, :]))
end
println("simulation complete")

#ANIMATION
    x = range(0, stop=L, length=nx)
    y = range(0, stop=L, length=nx)
    z = range(0, stop=L, length=nx)
    # Create all combinations of x, y, z coordinates
    x_coords = vec([xi for xi in x, yj in y, zk in z])
    y_coords = vec([yj for xi in x, yj in y, zk in z])
    z_coords = vec([zk for xi in x, yj in y, zk in z])

    anim = @animate for t in 1:speed:nt  
    tsec = (t-1)/fps
    scatter3d(x_coords, y_coords, z_coords, 
        marker_z=U[t, :, :, :][:],
        markersize=2,
        markerstrokewidth=0,
        c=:thermal,
        title="3D heat diffusion — t=$(round(tsec, digits=1)) s",
        xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",
        xlims=(0,L), ylims=(0,L), zlims=(0,L),
        colorbar=true, colorbar_title="Temperature [K]",
        clims=(minimum(U), maximum(U)),
        aspect_ratio=:equal,
        size=(600,600))
end
mp4(anim, "ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT\\SOLN_3D_Conduction.mp4", fps = fps)

#PLOT AVERAGE TEMPERATURE OVER TIME
time = range(0, stop=t_final, length=length(avg_temp))
p = plot(time, avg_temp,
    xlabel="Time [s]", ylabel="Average Temperature [K]",
    title="Average Temperature vs Time",
    linewidth=2, legend=false,
    size=(600,400))
savefig(p, "ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT\\SOLN_3D_Conduction.png")