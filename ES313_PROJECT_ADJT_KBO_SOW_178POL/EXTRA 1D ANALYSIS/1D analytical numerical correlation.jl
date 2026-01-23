    try
        using Plots
        using Statistics
        using Distributions
        using Printf
    catch
        import Pkg
        Pkg.add("Plots"); using Plots
        Pkg.add("Statistics"); using Statistics
        Pkg.add("Distributions"); using Distributions
        Pkg.add("Printf"); using Printf
    end


    # PARAMETERS
    
    α = 111e-6          # thermal diffusivity [m²/s]
    L = 1.0             # length of the rod [m]
    Ti = 298.15         # initial temperature [K]
    Tb = 1000.0         # boundary temperature [K]
    t_final = 600.0     # simulation time [s]
    nx = 65             # number of nodes
    fps = 10            # frames per second
    speed = 10          # plotting interval
    

    
    # ANALYTICAL SOLUTION

        Tavg(t) = Tb + (8 * (Ti - Tb) / π^2) * exp(-α * (π / L)^2 * t)
        t_analytical = range(0, t_final, length=500)
        T_analytical = Tavg.(t_analytical)


    # NUMERICAL SOLUTION

        dx = L / (nx - 1)
        dt = 1 / fps
        β = α * dt / dx^2
        if β >= 0.5
            error("β = $β too large for stability; increase nx or decrease fps")
        end
        nt = Int(t_final * fps)
        U = zeros(Float64, nt, nx)

        U[1, :] .= Ti
        U[:, 1] .= Tb
        U[:, end] .= Tb

        avg_temps = Float64[]
        time_points = Float64[]

        for t in 1:(nt-1)
            for i in 2:(nx-1)
                U[t+1, i] = U[t, i] + β * (U[t, i+1] - 2*U[t, i] + U[t, i-1])
            end
            
            if t % speed == 0
                push!(avg_temps, mean(U[t, :]))
                push!(time_points, (t - 1) * dt)
            end
        end


    # PLOT

    plot(t_analytical, T_analytical,
        label = "Analytical",
        xlabel = "Time [s]",
        ylabel = "Average Temperature [K]",
        title = "Average Temperature vs Time: Analytical vs Numerical",
        linewidth = 2,
        legend = :bottomright,
        grid = true,
        ylims = (Ti, Tb),   size = (800, 500))

    plot!(time_points, avg_temps,
        label = "Numerical",
        linewidth = 2,
        linestyle = :dash)

    
    annotate!([(t_final * 0.65, Ti + (Tb - Ti) * 0.85,
        text("Parameters:\n" *
                "Tᵢ = $(Float64(Ti)) K\n" *
                "Tբ = $(Float64(Tb)) K\n" *
                "L = $L m\n" *
                "dx = $(round(dx, digits=4)) m\n" *
                "dt = $(round(dt, digits=2)) s\n" *
                "β = $(round(β, digits=3))",
                :left, 9))])

    # STATISTICAL ANALYSIS

    t_start_analysis = 0.0
    t_end_analysis = t_final
    t_start_analysis = max(0.0, min(t_start_analysis, t_final))
    t_end_analysis = max(t_start_analysis, min(t_end_analysis, t_final))

    # Filter analytical data for the selected interval
    mask_analytical = (t_analytical .>= t_start_analysis) .& (t_analytical .<= t_end_analysis)
    t_analytical_filtered = t_analytical[mask_analytical]
    T_analytical_filtered = T_analytical[mask_analytical]

    numerical_interp = zeros(length(t_analytical_filtered))
    for (i, t) in enumerate(t_analytical_filtered)
        idx = searchsortedfirst(time_points, t)
        if idx > length(time_points)
            numerical_interp[i] = avg_temps[end]
        elseif idx == 1
            numerical_interp[i] = avg_temps[1]
        else
            # Linear interpolation
            t1, t2 = time_points[idx-1], time_points[idx]
            T1, T2 = avg_temps[idx-1], avg_temps[idx]
            numerical_interp[i] = T1 + (T2 - T1) * (t - t1) / (t2 - t1)
        end
    end

    errors = numerical_interp .- T_analytical_filtered
    abs_errors = abs.(errors)

    mae = mean(abs_errors)
    rmse = sqrt(mean(errors.^2))
    mean_err = mean(errors)
    std_err = std(errors)
    corr = cor(numerical_interp, T_analytical_filtered)

    n = length(errors)
    t_stat = mean_err / (std_err / sqrt(n))
    df = n - 1
    p_value = 2 * (1 - cdf(TDist(df), abs(t_stat)))

    open("ES313_PROJECT_ADJT_KBO_SOW_178POL\\EXTRA 1D ANALYSIS/1D analytical numerical comparison.txt", "w") do io
        println(io, "="^60)
        println(io, "STATISTICAL ANALYSIS")
        println(io, "="^60)
        println(io, "Mean Error:           ", round(mean_err, digits=4), " K")
        println(io, "Std Dev of Error:     ", round(std_err, digits=4), " K")
        println(io, "MAE:                  ", round(mae, digits=4), " K")
        println(io, "RMSE:                 ", round(rmse, digits=4), " K")
        println(io, "Correlation:          ", round(corr, digits=6))
        println(io, "="^60)
        println(io, "PAIRED T-TEST")
        println(io, "="^60)
        println(io, "H₀: Mean difference = 0")
        println(io, "H₁: Mean difference ≠ 0")
        println(io, "Time interval:        ", round(t_start_analysis, digits=2), " - ", round(t_end_analysis, digits=2), " s")
        println(io, "Degrees of freedom:   ", df)   
        println(io, "t-statistic:          ", round(t_stat, digits=4))
        println(io, "p-value:              ", @sprintf("%.4e", p_value))    
        println(io, "Result (α=0.05):      ", p_value > 0.05 ? "Fail to reject H₀" : "Reject H₀")
        println(io, "="^60)
    end

    savefig("ES313_PROJECT_ADJT_KBO_SOW_178POL\\EXTRA 1D ANALYSIS/1D analytical numerical comparison.png")

