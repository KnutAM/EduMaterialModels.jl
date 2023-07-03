### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 48db4cd0-15f0-11ee-0796-d3fffb0d18a2
begin
	io_log = open(joinpath(@__DIR__, "pkg_build.log"), "w")
	std_err_old = stderr
	std_out_old = stdout
	redirect_stdio(stdout=io_log, stderr=io_log)
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
		Pkg.PackageSpec(;name="PlutoUI"),
		Pkg.PackageSpec(;name="CairoMakie"),
		Pkg.PackageSpec(;name="Tensors"),
		Pkg.PackageSpec(;url="https://github.com/KnutAM/Newton.jl.git"),
		Pkg.PackageSpec(;url="https://github.com/KnutAM/MaterialModelsBase.jl.git"),
		Pkg.PackageSpec(;url="https://github.com/KnutAM/EduMaterialModels.jl.git"),
	])
	using PlutoUI, Tensors, MaterialModelsBase, EduMaterialModels
	import CairoMakie as CM
	redirect_stdio(stdout=std_out_old, stderr=std_err_old)
	close(io_log)
end;

# ╔═╡ 4b6d68c1-2b9f-45fc-bfda-a525cdf8679f
md"**Note:** If you are viewing this as a static html-page and want to change parameters: you must either open the notebook in [binder](https://binder.plutojl.org/v0.19.12/open?url=https%253A%252F%252Fknutam.github.io%252FEduMaterialModels.jl%252Fdev%252Fpluto_notebooks%252Fplasticity.jl) (very long loading time) or [download the file](https://knutam.github.io/EduMaterialModels.jl/dev/pluto_notebooks/plasticity.jl) and open it using [Pluto.jl](https://plutojl.org/)."

# ╔═╡ 2a7ab774-12b0-4a75-8266-891d8411c1f2
function plot_response(;Δϵ=0.01, num_steps=100, E=200.e3, Y0=200.0, Hiso=10e3, κ∞=100.0, Hkin=30.e3, β∞=100.0, use_explicit=false)
	tf(x) = SymmetricTensor{2,1}(tuple(x))
	ϵ_a = [tf(x) for x in range(0,Δϵ,num_steps)]
	ϵ_b = [tf(x) for x in range(Δϵ, -Δϵ, 2*num_steps)[2:end]]
	ϵ_c = [tf(x) for x in range(-Δϵ,0,num_steps)[2:end]]
	ϵ = append!(ϵ_a, ϵ_b, ϵ_c)
	m = EduMaterialModels.J2Plasticity(;
		e=EduMaterialModels.LinearIsotropicElasticity(;E=E, ν=0.3), 
		Y0=Y0, Hiso=Hiso, κ∞=κ∞, Hkin=Hkin, β∞=β∞)
	mw = use_explicit ? EduMaterialModels.ExplicitPlasticity(m) : m
	stress_state = UniaxialStress()
	σ = EduMaterialModels.simulate_response(mw, stress_state, ϵ, range(0,1,length(ϵ)))
	fig = CM.Figure()
	ax = CM.Axis(fig[1,1]; xlabel="ϵ₁₁ [%]", ylabel="σ₁₁ [MPa]")
	CM.lines!(ax, 100*first.(ϵ), first.(σ))
	CM.xlims!(ax, -100*Δϵ, 100*Δϵ)
	CM.ylims!(ax, -1000, 1000)
	return fig
end;

# ╔═╡ 7dc5f83b-7a07-4741-9a2e-42e9246a8999
begin
	E_slider = @bind E Slider(50:50:450.; show_value=true, default=200.)
	Y0_slider = @bind Y0 Slider(100.:50:500; show_value=true, default=250.)
	Hi_slider = @bind Hiso Slider(vcat(0, [2^n for n in 0:8]); show_value=true, default=32)
	κ∞_slider = @bind κ∞ Slider(vcat(50.:50:500, Inf); show_value=true, default=200.)
	Hk_slider = @bind Hkin Slider(vcat(0, [2^n for n in 0:8]); show_value=true, default=32)
	β∞_slider = @bind β∞ Slider(vcat(50.:50:500, Inf); show_value=true, default=200.)
	N_slider = @bind num_steps Slider([2^n for n in 3:10]; show_value=true, default=64)
	method_cb = @bind method CheckBox()
	md"""
	# Plasticity modeling
	## Intended learning outcomes
	1. Understand how the material parameters in the standard rate-independent Chaboche model affect the cyclic response of a material. Specifically, how kinematic hardening and isotropic hardening differ.
	2. How do the time discretization and time integration affect the material response, i.e. the number of time steps and explicit versus implicit algorithms. 

	## Model theory
	The plasticity model investigated in this notebook assumes linear elasticity via the additive decomposition, ``\boldsymbol{\epsilon} = \boldsymbol{\epsilon}_\mathrm{p} + \boldsymbol{\epsilon}_\mathrm{e}``, where the elastic strains, ``\boldsymbol{\epsilon}_\mathrm{e}`` are giving the stress.
	
	``
	\boldsymbol{\sigma} = \mathsf{E}:[\boldsymbol{\epsilon} - \boldsymbol{\epsilon}_\mathrm{p}], \quad \mathsf{E}=2G \mathsf{I}^\mathrm{dev} + 3K \boldsymbol{I}\otimes\boldsymbol{I}
	``
	
	where the elastic stiffness is assumed isotropic and is given by the shear and bulk moduli, ``G`` and ``K``. Young's modulus is given as ``E=9KG/(3K+G)``. 
	
	The yield criterion is based on the von Mises effective stress, 
	
	``
	\varPhi = f_\mathrm{vM}(\boldsymbol{\sigma}-\boldsymbol{\beta}) - [Y_0 + \kappa] , \quad
	f_\mathrm{vM}(\boldsymbol{x}) = \sqrt{\frac{3}{2}\boldsymbol{x}^\mathrm{dev}:\boldsymbol{x}^\mathrm{dev}}
	``
	
	where ``Y_0`` is the initial yield limit. The evolution law for the plastic strains is 
	
	``
	\dot{\boldsymbol{\epsilon}}_\mathrm{p} = \dot{\lambda} \boldsymbol{\nu}, \quad \boldsymbol{\nu} = \frac{\partial \varPhi}{\partial \boldsymbol{\sigma}}=\frac{3}{2} \frac{\boldsymbol{\sigma}^\mathrm{dev}-\boldsymbol{\beta}^\mathrm{dev}}{f_\mathrm{vM}(\boldsymbol{\sigma}-\boldsymbol{\beta})}
	``
	
	The hardening is governed by the Voce (isotropic) 
	
	``
	\dot{\kappa} = H_\mathrm{iso} \dot{\lambda} \left[1 - \frac{\kappa}{\kappa_\infty}\right]
	``

	and the Armstrong-Frederick (kinematic),
	
	``
	\dot{\boldsymbol{\beta}} = H_\mathrm{kin} \dot{\lambda} \frac{2}{3} \left[\boldsymbol{\nu} - \frac{3}{2}\frac{\boldsymbol{\beta}}{\kappa_\infty}\right]
	``

	laws, with isotropic and kinematic hardening modulii, ``H_\mathrm{iso}`` and ``H_\mathrm{kin}``, and isotropic and kinematic hardening saturation stresses, ``\kappa_\infty`` and ``\beta_\infty``. Finally, the KKT-conditions are used to obtain a rate-independent response,
	
	``
	\varPhi \leq 0, \quad \dot{\lambda} \geq 0, \quad \varPhi \dot{\lambda} = 0
	``´
	
	## Simulation parameters
	This notebook simulates the uniaxial stress response for ``\epsilon_{11}=\pm 1 \%`` during one cycle. You can adjust the material parameters, ``E``, ``Y_0``, ``H_\mathrm{iso}``, ``\kappa_\infty``, ``H_\mathrm{kin}``, and ``\beta_\infty``. In addition, you can choose how many steps are taken for each quarter cycle, as well as if explicit (Forward Euler) time integration should be used instead of implicit (Backward Euler). 
	
	| param | value        | unit |
	| ----- | ------------ | ---- |
	| E     | $(E_slider)  | GPa  |
	| Y0    | $(Y0_slider) | MPa  |
	| Hiso  | $(Hi_slider) | GPa  |
	| κ∞    | $(κ∞_slider) | MPa  |
	| Hkin  | $(Hk_slider) | GPa  |
	| β∞    | $(β∞_slider) | MPa  |
	| steps | $(N_slider)  | -    |  
	| explicit | $(method_cb) |   |
	
	"""
end

# ╔═╡ c5d0dd04-241e-48a8-944d-dfa61edeface
plot_response(;E=E*1e3, Y0=Y0, Hiso=1e3*Hiso, κ∞=κ∞, Hkin=1e3*Hkin, β∞=β∞, num_steps=num_steps, use_explicit=method)

# ╔═╡ Cell order:
# ╟─4b6d68c1-2b9f-45fc-bfda-a525cdf8679f
# ╟─48db4cd0-15f0-11ee-0796-d3fffb0d18a2
# ╟─2a7ab774-12b0-4a75-8266-891d8411c1f2
# ╟─7dc5f83b-7a07-4741-9a2e-42e9246a8999
# ╟─c5d0dd04-241e-48a8-944d-dfa61edeface
