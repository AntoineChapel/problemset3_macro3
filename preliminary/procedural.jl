using Plots
α = 0.40
β = 0.60
θ = 10
k_ss = (1/(α*β*θ))^(1/(α-1))
c_ss = 10*(k_ss^α) - k_ss
#transversality condition:
t = 10 .^collect(0:5)
k = 10 .^collect(0:5)
value = β.^t * transpose(log.(θ*(k.^α)) .+ log(1-α*β))

# => Whatever the value of k, the transversality condition is satisfied when t goes to infinite

g_ana(k) = α*β*θ*(k^α)
F(k) = θ*(k^α)
C(k) = F(k) - g_ana(k)
nk = 30
kmin, kmax = 2, 15
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
euler_adj = euler[2:end]
println(intertemporal)
println(euler_adj)


#both arrays are approximately equal => The Euler equation is verified with the policy function




#c)

k_ss_analytical = (1/(α*β*θ))^(1/(α-1))

nk = 5
α = 0.4
β = 0.6
θ = 10
itermax = 10
kmin, kmax = 2, 10
kg = collect(kmin:(kmax-kmin)/(nk-1):kmax)
u(ct) = log(ct)
F(kt) = θ*(kt^α)


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

plot(kt, kp_opt[:,end], label="Iterated", lw=2)
#It seems that the agent saves a fixed amount for some values of k, but this fixed value increases as capital accumulates. We would obtain a more precise result with more iterations and a larger grid.

plot!(kt, g_ana.(kt), label="Analytical", lw=2)
#the result of the iterative algorithm is an approximation of the analytical version

plot(kt, CF, label="Consumption", lw=2)
#Consumption increases when capital accumulates, which makes sense especially since there is no depreciation of capital

k_opti = kp_opt[:, end]
diff = kp_opti - kt
kss_simu = kt[findmin(diff.^2)[2]]
println("The difference between the fixed point obtained analytically and through the iterative procedure is: $(kss_simu-k_ss_analytical)")




#e)

tol = 1E-5
iterinit = 1
normdiff = Inf
maxiter = 1000


V0 = ones(nk, 1)
Vi = Array{Float64}(undef, nk, maxiter)
V = hcat(V0, Vi)
PF = Array{Float64}(undef, nk, maxiter)
nk = 5
Objgrid = Array{Float64}(undef, nk, nk, maxiter)
kmin, kmax = 2, 15
kg = collect(kmin:((kmax-kmin)/(nk-1)):kmax)


u(ct) = log(ct)
F(kt) = θ*(kt^α)

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

V_0 = V0
V_1 = Array{Float64}(undef, nk, 1)
iteration = iterinit

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
#plot!(kt, kt, legend=false, lw=2)

kp_opti = kp_opt[:, iteration-1]
diff = kp_opti - kt
kss_simu = kt[findmin(diff.^2)[2]]

println("The difference between the fixed point obtained analytically and through the iterative procedure is: $(kss_simu-k_ss_analytical)")
#Clearly, the result obtained is better through this procedure, for an identical number of grid elements

plot(kt, kp_opt[:, iteration-1], label="Final Policy function", lw=2)

value_function_final = V[:, iteration]
plot!(kt, value_function_final, label="Final Value function", lw=2)


#f)

kg = collect(0.05:0.05:10)
nk = 200
tol = 1E-5
iterinit = 1
normdiff = Inf
maxiter = 1000


V0 = ones(nk, 1)
Vi = Array{Float64}(undef, nk, maxiter)
V = hcat(V0, Vi)
PF = Array{Float64}(undef, nk, maxiter)
Objgrid = Array{Float64}(undef, nk, nk, maxiter)


u(ct) = log(ct)
F(kt) = θ*(kt^α)

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

V_0 = V0
V_1 = Array{Float64}(undef, nk, 1)
iteration = iterinit

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

plot(kt, kp_opt[:, iteration-1], label="Optimal Policy function", lw=2)
plot!(kt, g_ana.(kt), label="Analytical Policy function", lw=1.5)
#plot!(kt, kt, legend=false, lw=2)

kp_opti = kp_opt[:, iteration-1]
diff = kp_opti - kt
kss_simu = kt[findmin(diff.^2)[2]]

println("The difference between the fixed point obtained analytically and through the iterative procedure is: $(kss_simu-k_ss_analytical)")

#the error term is even smaller because the grid is thinner
