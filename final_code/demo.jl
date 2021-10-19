include("myenv\\src\\Exo5_functions.jl")


n=100

grid = grid_generator(2, 15, n)
k_simu, k_analytical, kp_simu, kp_analytical, CF = valuefunctioniteration_itermax(0.4, 0.6, 10, grid, 30)
plot_k_results(grid, kp_simu, kp_analytical)
plot_c_results(grid, CF)


grid = grid_generator(2, 15, n)
kss_ana_conv, kss_simu_conv, k_ana_conv, k_simu_conv, CF_conv, V_conv, v_final = valuefunctioniteration_convergence(0.4, 0.6, 10, grid, n, "square", 1E-5, 10)



plot_k_results(grid, k_simu_conv, k_ana_conv)
plot_c_results(grid, CF_conv)
plot(grid, v_final, label="Final Value function", lw=2)
