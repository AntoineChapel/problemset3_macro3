using Plots

α = 0.40
β = 0.60
θ = 10


k_ss = (1/(α*β*θ))^(1/(α-1))
c_ss = 10*(k_ss^α) - k_ss


#b) transversality condition, Euler equations:

g_ana(k) = α*β*θ*(k^α)

function generate_path(k0, g, n)
    k_hist = zeros(n)
    k_hist[1] = k0
    for i in 2:n
        k_hist[i] = g(k_hist[i-1])
    end
    return k_hist
end


k_hist = generate_path.(10 .^-5:15, g_ana, 10)
plot(k_hist, legend=false)


nk = 10
kmin, kmax = 2, 15
k_ss = (1/(α*β*θ))^(1/(α-1))
c_ss = θ*(k_ss^α) - k_ss
g_ana(k) = α*β*θ*(k^α)
F(k) = θ*(k^α)
C(k) = F(k) - g_ana(k)
kg = collect(kmin:((kmax-kmin)/(nk-1)):kmax)
c = zeros(nk)
intertemporal = zeros(nk-1)
kp = zeros(nk)
k_hist = zeros(nk)
k_hist[1] = 1


for i in 2:nk
    k_hist[i] = g_ana(k_hist[i-1])
end
c_hist = C.(k_hist)

for j in 1:nk-1
    intertemporal[j] = c_hist[j+1]/c_hist[j]
end
euler = α*β*θ*k_hist.^(α - 1)
euleradjusted = euler[2:end]
println(intertemporal)
println(euleradjusted)


#both arrays are approximately equal => The Euler equation is verified with the policy function



#c), d)

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

    CF = zeros(nk, itermax)
    kt = copy(kg)
    kp_opt = zeros(nk, itermax)

    for iteration in 1:itermax
        for j in 1:nk
            kp_opt[j, iteration] = kg[Int(PF[j, iteration])]
            CF[j, iteration] = F(kt[j]) - kp_opt[j, end]
        end
    end

    kp_opti = kp_opt[:, end]
    diff = kp_opti - kt
    kss_simu = kt[findmin(diff.^2)[2]]
    kp_opt_analytical = pf_analytical.(grid)
    println("The difference between the fixed point obtained analytically and through the iterative procedure is: $(kss_simu-k_ss_analytical)")
    return kss_simu, k_ss_analytical, kp_opt, kp_opt_analytical, CF, V
end

n=5
grid = grid_generator(2, 10, n)
k_simu, k_analytical, kp_simu, kp_analytical, CF, V = valuefunctioniteration_itermax(0.4, 0.6, 10, grid, 30)

plot(grid, kp_simu[:, 2:end], legend=false , title="Policy Function iteration", lw=2)
plot(grid, CF, legend=false, title = "Consumption Function iteration", lw=2)
plot_k_results(grid, kp_simu, kp_analytical)
plot(grid, V[:, 2:end],legend=false, title="Value Function iteration", lw=1)


v_analytical(α, β, θ, k) = ((1/(1 - β)) * (log(θ - α*β*θ) + α*β/(1-α*β)*log(α*β*θ)))  + ((α/(1-α*β))*log(k))
v_ss = v_analytical.(0.4, 0.6, 10, grid)
plot!(grid, v_ss, label="analytical value function", lw=5, color="black")




#e)
#press ctrl+j

using Plots
k_ss = (1/(α*β*θ))^(1/(α-1))

tol = 1E-5
iterinit = 1
normdiff = Inf
maxiter = 1000


#for e):

#nk = 5
#kmin, kmax = 2, 10

#for f):
nk = 200
kmin, kmax = 0.05, 10



V0 = ones(nk, 1)
Vi = Array{Float64}(undef, nk, maxiter)
V = hcat(V0, Vi)
PF = Array{Float64}(undef, nk, maxiter)
Objgrid = Array{Float64}(undef, nk, nk, maxiter)
kg = collect(kmin:((kmax-kmin)/(nk-1)):kmax)


u(ct) = log(ct)
F(kt) = θ*(kt^α)
g_ana(k) = α*β*θ*(k^α)


V_0 = V0
V_1 = Array{Float64}(undef, nk, 1)
iteration = iterinit
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
while normdiff >= tol && iteration <= maxiter
    for i in 1:nk
        for j in 1:nk
            Objgrid[i, j, iteration] = Ut[i, j] + β*V[j, iteration]
        end
        V[i, iteration+1], PF[i, iteration] = findmax(Objgrid[i, :, iteration])
    end

    global V_1 = V[:, iteration+1]
    global normdiff = maximum(abs.(V_1-V_0))
    global V_0 = copy(V_1)
    iteration += 1
    println(iteration)
end


kt = copy(kg)
CF = Array{Float64}(undef, nk)
kp_opt = Array{Float64}(undef, nk, maxiter)

for iter in 1:iteration-1
    for j in 1:nk
        kp_opt[j, iter] = kg[Int(PF[j, iter])]
        CF[j] = F(kt[j]) - kp_opt[j, iter]
    end
end

plot(kt, kp_opt[:, 2:iteration-1], legend=false, lw=2)

kp_opti = kp_opt[:, iteration-1]
diff = kp_opti - kt
kss_simu = kt[findmin(diff.^2)[2]]

println("The difference between the fixed point obtained analytically and through the iterative procedure is: $(kss_simu-k_ss)")
#Clearly, the result obtained is better through this procedure, for an identical number of grid elements

plot(kt, kp_opt[:, iteration-1], label="Final Policy function", legend=:bottomright, lw=2)
plot!(kt, g_ana.(kt), label="Analytical policy function")


value_function_final = V[:, iteration]
plot(kt, value_function_final, label="Final iterated Value function", lw=2, legend=:bottomright, color="red")
plot!(kt, v_analytical.(0.4, 0.6, 10, kt), label="Analytical value function", lw=1, color="black")