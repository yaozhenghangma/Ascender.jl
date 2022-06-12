#%% an example from https://github.com/tamuhey/python_1d_dft/
using Plots
using LinearAlgebra
using Formatting
gr()

#%% grid
number_grid = 200
x = range(-5, 5, length=number_grid)
y = sin.(x)
plot(x, y)


#%% first order differentiation
h = x[2] - x[1]
D = diagm(1=>ones(number_grid-1)) - I
D = D./h;

#%% second order differentiation
D2 = diagm(1=>ones(number_grid-1)) + diagm(-1=>ones(number_grid-1)) - 2*I
D2 = D2./h^2;

#%%
plot(x, y, lable="f")
plot!(x[1:end-1], (D*y)[1:end-1], label="df")
plot!(x[2:end-1], (D2*y)[2:end-1], label="d2f")

#%% non-interacting electrons
eigenvalue, eigenvector1 = eigen(-D2./2)
plot()
for i in 1:5
    plot!(x, eigenvector1[:,i], label=format(eigenvalue[i], precision=4))
end
plot!()

#%% harmonic oscillator
X = diagm(0 => x.*x)
eigenvalue, eigenvector2 = eigen(-D2./2 + X)
plot()
for i in 1:5
    plot!(x, eigenvector2[:,i], label=format(eigenvalue[i], precision=4))
end
plot!()

#%% well potential
w = fill(1e10, size(x))
@. w[2>x>-2] = 0.0
w = diagm(0 => w)
eigenvalue, eigenvector3 = eigen(-D2./2 + w)
plot()
for i in 1:5
    plot!(x, eigenvector3[:,i], label=format(eigenvalue[i], precision=4))
end
plot!()


#%% density
function integral(x, y; dims=1)
    dx = x[2]-x[1]
    return sum(y.*dx, dims=dims)
end
number_electrons = 17
function get_nx(number_electrons, ψ, x)
    factor = integral(x, abs2.(ψ), dims=1)
    normalized_ψ = ψ./sqrt.(factor)

    fn = [2 for _ in 1:number_electrons÷2]
    if (number_electrons%2) == 1
        push!(fn, 1)
    end

    res = zeros(size(ψ)[1])
    for i in zip(fn, 1:size(normalized_ψ)[2])
        res += i[1].*abs2.(normalized_ψ[:, i[2]])
    end
    return res
end

plot(get_nx(number_electrons, eigenvector1, x), label="non")
plot!(get_nx(number_electrons, eigenvector2, x), label="harm")
plot!(get_nx(number_electrons, eigenvector3, x), label="well")


#%% exchange energy
function get_exchange(nx, x)
    energy = -3/4 .* (3/π)^(1/3) .* integral(x, nx.^(4/3))
    potential = -(3/π)^(1/3) .* nx .* (1/3)
    return energy, potential
end

#%% Coulomb potential
function get_hatree(nx, x; eps=1e-1)
    h = x[2]-x[1]
    denominator_factor = 1.0./
        sqrt.((x*transpose(ones(size(x))) - ones(size(x))*transpose(x)).^2 .+ eps)
    energy = sum(nx*transpose(nx).*0.5.*h^2 .* denominator_factor, dims=[1,2])
    potential = dropdims(
        sum(ones(size(nx))*transpose(nx) .* h .* denominator_factor, dims=2),
        dims=2)
    return energy, potential
end

#%% self-consistency calculation
function print_log(i, log)
    printfmtln("step: {:<5d} energy: {:<10.4f} energy_diff: {:<f}",
        i, log["energy"][i], log["energy_diff"][i])
    return nothing
end

max_iter = 100
energy_tolerance = 1e-5
log = Dict("energy"=>Array{Float64, 1}([]), "energy_diff"=>Array{Float64, 1}([]))

nx = zeros(number_grid)
energy = zeros(size(x))
ψ = zeros(size(eigenvector1))
for i in 1:max_iter
    ex_energy, ex_potential = get_exchange(nx, x)
    ha_energy, ha_potential = get_hatree(nx, x, eps=0.1)
    #print(nx[1:3])

    H = -D2./2 + diagm(0=> ex_potential+ha_potential + x.*x)

    global energy, ψ = eigen(H)

    push!(log["energy"], energy[1])
    if isempty(log["energy_diff"])
        energy_diff = Inf
    else
        energy_diff = energy[1] - log["energy"][end-1]
    end
    push!(log["energy_diff"], energy_diff)
    print_log(i, log)

    if abs(energy_diff) < energy_tolerance
        println("converged!")
        break
    elseif i == max_iter
        println("not converged!")
    end
    global nx = get_nx(number_electrons, ψ, x)
end

plot()
for i in 1:5
    plot!(x, ψ[:,i], label=format(energy[i], precision=4))
end
plot!()

#%%
plot(nx)
plot!(get_nx(number_electrons, eigenvector2, x), label="non")
