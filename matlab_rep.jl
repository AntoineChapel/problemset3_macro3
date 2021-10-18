using LinearAlgebra, Statistics, Plots


α = 0.4
β = 0.6
δ = 1
θ = 10


nk = 50
itermax = 10
u(ct) = log(ct)
F(kt) = θ*(kt^α)

k_ss_analytical = ((θ*α)/((1/β)-(1-δ)))^(1/(1-α))
y_ss = F(k_ss_analytical)
c_ss = y_ss - δ*k_ss_analytical

kmin, kmax = 2, 15
kg = collect(kmin:((kmax-kmin)/(nk-1)):kmax)


Ut = Array{Float64}(undef, nk, nk)
for i in 1:nk
    for j in 1:nk
        k = kg[i]
        kp = kg[j]
        c = F(k) + (1 - δ)*k - kp

        if kp < 0 || c < 0
            Ut[i, j] = -Inf
        else
            Ut[i, j] = u(c)
        end
    end
end


V0 = ones(nk, 1)
Vi = zeros(nk, itermax)
V = hcat(V0, Vi)
PF = zeros(nk, itermax)

Objgrid = zeros(nk, nk, itermax)
#Objgrid = zeros(nk, nk, 0)


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
        CF[j] = F(kt[j]) + (1 - δ)*kt[j] - kp_opt[j, end]
    end
end

plot(kt, kp_opt[:, 2:end], legend=false, lw=2)
plot!(kt, kt, legend=false, lw=2)

kp_opti = kp_opt[:, end]

diff = kp_opti - kt
kss_simu = kt[findmin(diff.^2)[2]]
k_ss_analytical

kss_simu
k_ss_analytical

kp_opt