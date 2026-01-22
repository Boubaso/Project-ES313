try
    using Plots
    using Statistics
catch
    import Pkg; Pkg.add("Plots"); using Plots
    import Pkg; Pkg.add("Statistics"); using Statistics
end
# INPUT
α = 111*10^(-6)     # thermal diffusivity m^2/s
L = 1.0             # length of the surf [m]
nx = 65             # number of knodes
t_final = 600
speed = 10       # speed factor for simulation
# CODE
fps = 10
dx = L/(nx-1)
dt = 1/fps
β = α * dt / (dx^2)
if β >= 0.05 error("β too large for stability condition; increase nx or decrease fps") end
nt = t_final*fps
U_vec = []
U = ones(Float64, nt, nx, nx)*298.15
U[:, 1, :] .= 1000
U[:, end, :] .= 1000
function heat_diffusion_2D(U, β, nx)
     U_new = copy(U)
     for i in 2:(nx-1)
         for j in 2:(nx-1)
             U_new[i, j] = U[i, j] + β * (U[i+1, j] + U[i-1, j] - 4*U[i, j] + U[i, j+1] + U[i, j-1])
         end
     end
     return U_new
end
for t in 1:(nt-1)
    U[t+1, :, :] = heat_diffusion_2D(U[t, :, :], β, nx)
end
clims = (minimum(U), maximum(U))
x = range(0, L, length=nx)
y = x
anim = @animate for k in 1:speed:nt
    tsec = (k-1)/fps
    heatmap(x, y, U[k, :, :]', clims=clims, c=:thermal,
            title="2D heat diffusion — t=$(round(tsec, digits=1)) s",
            xlabel="x [m]", ylabel="y [m]", aspect_ratio=1, size=(600,600),
            xlims=(0,L), ylims=(0,L),
            colorbar=true, colorbar_title="Temperature [K]")
end
mp4(anim, "ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT\\SOLN_2D_Conduction.mp4", fps = fps)

avg_temps = [mean(U[k, :, :]) for k in 1:nt]
time_vec = [(k-1)/fps for k in 1:nt]
plot(time_vec, avg_temps,
    xlabel="Time [s]", ylabel="Average Temperature [K]",
    title="Average Temperature vs Time",
    linewidth=2, legend=false,
    size=(600,400))
savefig("ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT\\SOLN_2D_Conduction.png")

