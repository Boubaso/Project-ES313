###### SETUP ######
try
    using Plots
catch
    import Pkg; Pkg.add("Plots")
    using Plots
end

println("Packages loaded successfully.")
println("Starting FK parameter calculation...")
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
# other 
    Tmin = 350.0    # K
    Tmax = 800.0    # K          

###### FUNCTIONS ######

# Frank–Kamenetskii parameter
    function FK_parameter(T, L, material)
        return (material.Q * material.A * L^2) /(material.ρ * material.cp * material.α * T^2) * exp(-material.E / (material.Rg * T))
    end
# Dakholmer number
    function Dakholmer_number(geometry)
        if geometry == :plane
            δc = 0.8785
        elseif geometry == :cylinder
            δc = 1.0
        end
    end

# Critical temperature solver
    function critical_temperature(L, material, geometry; Tmin=350.0, Tmax=800.0, tol=1e-6, maxiter=100)
        δc = Dakholmer_number(geometry)
        f(T) = FK_parameter(T, L, material) - δc
        a, b = Tmin, Tmax
        fa, fb = f(a), f(b)
        if fa * fb > 0
            error("Root not bracketed — adjust Tmin/Tmax")
        end
        for _ in 1:maxiter
            c = 0.5 * (a + b)
            fc = f(c)
            abs(fc) < tol && return c
            if fa * fc < 0
                b, fb = c, fc
            else
                a, fa = c, fc
            end
        end
        error("Critical temperature did not converge")
    end

# Ignition time estimate
    function ignition_time(T, L, material)
        δ = FK_parameter(T, L, material)
        δ <= 1 && return Inf
        td = L^2 / material.α * (π^2 / 4) / (δ - 1)
        return td / sqrt(δ - 1)
    end

# Example usage
L = 1e-3   # 1 mm hot-spot radius
Tc_plane = critical_temperature(L, material, :plane)
Tc_cylinder = critical_temperature(L, material, :cylinder)
    println("Critical temperature Tc plane = $(round(Tc_plane, digits=1)) K")
    println("Critical temperature Tc cylinder = $(round(Tc_cylinder, digits=1)) K")
Ttest_plane = Tc_plane * 1.05
Ttest_cylinder = Tc_cylinder * 1.05
    println("Testing ignition at T = $(round(Ttest_plane, digits=1)) K (plane)")
    println("Testing ignition at T = $(round(Ttest_cylinder, digits=1)) K (cylinder)")
t_ign_plane = ignition_time(Ttest_plane, L, material)
t_ign_cylinder = ignition_time(Ttest_cylinder, L, material)
    println("Ignition time at 1.05 Tc (plane) ≈ $(round(t_ign_plane, digits=2)) s")
    println("Ignition time at 1.05 Tc (cylinder) ≈ $(round(t_ign_cylinder, digits=2)) s")

    # Save output to text file
    output_file = "ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT\\FK_analytical_output.txt"
    open(output_file, "w") do io
        println(io, "Frank-Kamenetskii Analysis Results")
        println(io, "================================================")
        println(io, "")
        println(io, "Material: RDX")
        println(io, "Hot-spot radius L = $(L*1e3) mm")
        println(io, "")
        println(io, "Critical temperature Tc plane = $(round(Tc_plane, digits=1)) K")
        println(io, "Critical temperature Tc cylinder = $(round(Tc_cylinder, digits=1)) K")
        println(io, "")
        println(io, "Testing ignition at T = $(round(Ttest_plane, digits=1)) K (plane)")
        println(io, "Testing ignition at T = $(round(Ttest_cylinder, digits=1)) K (cylinder)")
        println(io, "")
        println(io, "Ignition time at 1.05 Tc (plane) ≈ $(round(t_ign_plane, digits=2)) s")
        println(io, "Ignition time at 1.05 Tc (cylinder) ≈ $(round(t_ign_cylinder, digits=2)) s")
    end
    println("Results saved to $output_file")

# Plot f(T) = FK_parameter 
    T_range = range(Tmin, Tmax, length=1000)

    f_plane(T) = FK_parameter(T, L, material) - Dakholmer_number(:plane)
    f_cylinder(T) = FK_parameter(T, L, material) - Dakholmer_number(:cylinder)
    f_values_plane = [f_plane(T) for T in T_range]
    f_values_cylinder = [f_cylinder(T) for T in T_range]

    p_plane = plot(T_range, f_values_plane, 
        xlabel="Temperature T (K)", 
        ylabel="f(T) = δ - δ_c",
        label="f(T)",
        linewidth=2,
        title="Frank-Kamenetskii Parameter (Plane)",
        grid=true,
        legend=:topright)
    hline!([0], color=:red, linestyle=:dash, label="Critical threshold")
    vline!([Tc_plane], color=:green, linestyle=:dash, label="Tc = $(round(Tc_plane, digits=1)) K")


    p_cylinder = plot(T_range, f_values_cylinder, 
        xlabel="Temperature T (K)", 
        ylabel="f(T) = δ - $(Dakholmer_number(:cylinder))",
        label="f(T)",
        linewidth=2,
        title="Frank-Kamenetskii Parameter (Cylinder)",
        grid=true,
        legend=:topright)
    hline!([0], color=:red, linestyle=:dash, label="Critical threshold")
    vline!([Tc_cylinder], color=:green, linestyle=:dash, label="Tc = $(round(Tc_cylinder, digits=1)) K")
    
    display(plot(p_plane, p_cylinder, layout=(1,2), size=(1000,450)))

    # Save the plot
    savefig(p_plane, "ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT\\FK_analytical_plot_plane.png")
    println("Plot FK saved to ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT\\FK_analytical_plot_plane.png")
    savefig(p_cylinder, "ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT\\FK_analytical_plot_cylinder.png")
    println("Plot FK saved to ES313_PROJECT_ADJT_KBO_SOW_178POL\\OUTPUT\\FK_analytical_plot_cylinder.png")