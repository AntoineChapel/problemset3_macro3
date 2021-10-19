using Plots


function grid_generator(kmin, kmax, n)
    return collect(kmin:(kmax-kmin)/(n-1):kmax)
end


function plot_k_results(grid, k_simu, k_ana)
    plot(grid, k_simu[:,end], label="Iterated", lw=2)
    plot!(grid, k_ana, label="Analytical", lw=2)
    
end

function plot_c_results(grid, CF)
    plot(grid, CF, label="Consumption", lw=2)
end



function valuefunctioniteration_itermax(α, β, θ, kg, itermax=10)
    k_ss_analytical = (1/(α*β*θ))^(1/(α-1))
    c_ss = θ*(k_ss_analytical^α) - k_ss_analytical
    F(k) = θ*(k^α)
    C(k) = F(k) - pf_analytical(k)
    u(c) = log(c)
    pf_analytical(k) = α*β*θ*(k^α)
    nk = length(grid)
    Ut = Array{Float64}(undef, nk, nk)
    for i in 1:nk
        for j in 1:nk
            k = kg[i]
            kp = kg[j]
            c = F(k) - kp

            if kp < 0 || c < 0
                Ut[i, j] = -Inf
            else
                Ut[i, j] = u(c)
            end
        end
    end
    V0 = zeros(nk, 1)
    Vi = zeros(nk, itermax)
    V = hcat(V0, Vi)
    PF = zeros(nk, itermax)
    Objgrid = zeros(nk, nk, itermax)

    for iteration in 1:itermax
        for i in 1:nk
            for j in 1:nk
                Objgrid[i, j, iteration] = Ut[i, j] + β*V[j, iteration]
            end
            V[i, iteration+1], PF[i, iteration] = findmax(Objgrid[i, :, iteration])
        end
    end

    CF = zeros(nk)
    kt = copy(kg)
    kp_opt = zeros(nk, itermax)

    for iteration in 1:itermax
        for j in 1:nk
            kp_opt[j, iteration] = kg[Int(PF[j, iteration])]
            CF[j] = F(kt[j]) - kp_opt[j, end]
        end
    end

    
    kp_opti = kp_opt[:, end]
    diff = kp_opti - kt
    kss_simu = kt[findmin(diff.^2)[2]]
    kp_opt_analytical = pf_analytical.(grid)
    println("The difference between the fixed point obtained analytically and through the iterative procedure is: $(kss_simu-k_ss_analytical)")
    return kss_simu, k_ss_analytical, kp_opt, kp_opt_analytical, CF
end





function valuefunctioniteration_convergence(α, β, θ, grid, nk, metric="abs", tol=1E-5, maxiter=1000)
    
    iterinit = 1
    normdiff = Inf
    k_ss_analytical = (1/(α*β*θ))^(1/(α-1))

    V0 = ones(nk)
    Vi = Array{Float64}(undef, nk, maxiter)
    V = hcat(V0, Vi)

    kg = grid
    PF = Array{Float64}(undef, nk, maxiter)
    Objgrid = Array{Float64}(undef, nk, nk, maxiter)

    u(ct) = log(ct)
    F(kt) = θ*(kt^α)
    pf_analytical(k) = α*β*θ*(k^α)

    Ut = Array{Float64}(undef, nk, nk)
    for i in 1:nk
        for j in 1:nk
            k = kg[i]
            kp = kg[j]
            c = F(k) - kp
            if kp < 0 || c < 0
                Ut[i, j] = -Inf
            else
                Ut[i, j] = u(c)
            end
        end
    end

    iteration = iterinit

    while normdiff >= tol && iteration <= maxiter
        for i in 1:nk
            for j in 1:nk
                Objgrid[i, j, iteration] = Ut[i, j] + β*V[j, iteration]
            end
            V[i, iteration+1], PF[i, iteration] = findmax(Objgrid[i, :, iteration])
        end

        global V_n = V[:, iteration+1]
        global V_p = V[:, iteration]

        if metric == "abs"
            global normdiff = maximum(abs.(V_n-V_p))
        else
            global normdiff = maximum((V_n-V_p).^2)
        end
        iteration += 1
        println(iteration)
    end

    CF = Array{Float64}(undef, nk)
    kp_opt = Array{Float64}(undef, nk, maxiter)

    for iter in 1:iteration-1
        for j in 1:nk
            kp_opt[j, iter] = kg[Int(PF[j, iter])]
            CF[j] = F(grid[j]) - kp_opt[j, iter]
        end
    end

    value_function_final = V[:, iteration]

    kp_opti = kp_opt[:, iteration-1]
    diff = kp_opti - grid
    kss_simu = grid[findmin(diff.^2)[2]]

    k_analytical = pf_analytical.(grid)
    k_simu = copy(kp_opti)
    println("The difference between the fixed point obtained analytically and through the iterative procedure is: $(kss_simu-k_ss_analytical)")
    return k_ss_analytical, kss_simu, k_analytical, k_simu, CF, V, value_function_final

end





