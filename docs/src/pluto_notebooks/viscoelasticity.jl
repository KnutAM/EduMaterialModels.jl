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
		Pkg.PackageSpec(;url="https://github.com/KnutAM/EduMaterialModels.jl", rev="main"),
	])
	using PlutoUI, Tensors, MaterialModelsBase, EduMaterialModels
	import CairoMakie as CM
	redirect_stdio(stdout=std_out_old, stderr=std_err_old)
	close(io_log)
end;

# ╔═╡ 4b6d68c1-2b9f-45fc-bfda-a525cdf8679f
md"**Note:** If you are viewing this as a static html-page and want to change parameters: you must either open the notebook in [binder](https://binder.plutojl.org/v0.19.12/open?url=https%253A%252F%252Fknutam.github.io%252FEduMaterialModels.jl%252Fdev%252Fpluto_notebooks%252Fviscoelasticity.jl) (very long loading time) or [download the file](https://knutam.github.io/EduMaterialModels.jl/dev/pluto_notebooks/viscoelasticity.jl) and open it using [Pluto.jl](https://plutojl.org/)."

# ╔═╡ 2a7ab774-12b0-4a75-8266-891d8411c1f2
function plot_response(;t_ramp=1.0, E1=30e3, E2=1e3, η=1e-3, ϵ=0.001, num_steps=100)
	t_total=5.0
	tf(x) = SymmetricTensor{2,1}(tuple(x))
	ϵ_a = [tf(x) for x in range(0,ϵ,num_steps)]
	ϵ_b = [tf(ϵ) for _ in 2:num_steps]
	ϵ_vec = vcat(ϵ_a, ϵ_b)
	t = vcat(collect(range(0, t_ramp, num_steps)), range(t_ramp, t_total, num_steps)[2:end])
	
	m = EduMaterialModels.Zener1D(;E1, E2, η)
	stress_state = UniaxialStress()
	σ = EduMaterialModels.simulate_response(m, stress_state, ϵ_vec, t)
	
	fig = CM.Figure()
	ax = CM.Axis(fig[1,1]; xlabel="ϵ₁₁ [%]", ylabel="σ₁₁ [MPa]")
	CM.lines!(ax, 100*first.(ϵ_vec), first.(σ))
	CM.xlims!(ax, 0, 100*ϵ)
	CM.ylims!(ax, 0, 200)

	ax = CM.Axis(fig[2,1]; xlabel="t [s]", ylabel="σ₁₁ [MPa]")
	CM.lines!(ax, t, first.(σ))
	CM.xlims!(ax, 0,t_total)
	CM.ylims!(ax, 0, 200)
	
	return fig
end;

# ╔═╡ 7dc5f83b-7a07-4741-9a2e-42e9246a8999
begin
	E1_slider = @bind E1 Slider(0:5:100.; show_value=true, default=50.)
	E2_slider = @bind E2 Slider(0:5:100.; show_value=true, default=50.)
	η_slider = @bind η Slider(0:2:50.; show_value=true, default=2.0)
	tramp_slider = @bind t_ramp Slider([2.0^n for n in -4:2], show_value=true, default=1.0)
	N_slider = @bind num_steps Slider([2^n for n in 3:10]; show_value=true, default=64)
	
	md"""
	# Visco-elastic modeling
	## Intended learning outcomes
	1. Understand how the different material parameters of the Zener model interact with each other. 

	## Model theory
	The Zener model consists of a linear spring (stiffness ``E_1``), connected in parallel with a Maxwell element. The maxwell element consists of a linear spring (stiffness ``E_2``) connected in series with a damper (coefficient ``\eta``).

	![Zener model figure is missing](https://knutam.github.io/EduMaterialModels.jl/dev/pluto_notebooks/zener_rheology.svg)

	Based on this rheological model, we consider the total stress as the sum of each parallel connection
	
	``\sigma = \sigma_1 + \sigma_2``

	where

	``\sigma_1 = E_1 \epsilon`` 

	and 

	``\sigma_2 = E_2 \epsilon_\mathrm{e}``

	The elastic strain, ``\epsilon_\mathrm{e}``, stems from an additive decomposition of the strain,

	``\epsilon = \epsilon_\mathrm{e} + \epsilon_\mathrm{v}``

	where ``\epsilon_\mathrm{v}`` is the viscous strains. The consitutive relationship for the damper yields

	``\sigma_2 = \eta \dot{\epsilon}_\mathrm{v}``

	thus providing an evolution law for ``\epsilon_\mathrm{v}``,

	``\dot{\epsilon}_\mathrm{v} = \frac{E_2}{\eta} \epsilon_\mathrm{e}``

	(Noting that ``\epsilon_\mathrm{e}`` is a function of ``\epsilon_\mathrm{v}``. However, even for the Backward Euler implicit integration, it is fairly straightforward to find the analytical solution)

	## Simulation parameters
	This notebook simulates the uniaxial stress response for an initial ramp of ``\epsilon_{11}`` to ``0.1 \%`` strain, during the time ``t_\mathrm{ramp}``. You can adjust the material parameters, ``E_1``, ``E_2``, and ``\eta``. In addition, you can choose the ramping time and how many steps are taken (same number for  the ramp and hold time). 
	
	| param | value        | unit |
	| ----- | ------------ | ---- |
	| E1    | $(E1_slider) | GPa  |
	| E2    | $(E2_slider) | GPa  |
	| η     | $(η_slider)  | GPa/s|
	| ramp  | $(tramp_slider)  | s|
	| steps | $(N_slider)  | -    |  
	
	"""
end

# ╔═╡ c5d0dd04-241e-48a8-944d-dfa61edeface
plot_response(;E1=E1*1e3, E2=E2*1e3, η=η*1000, t_ramp, num_steps)

# ╔═╡ Cell order:
# ╟─4b6d68c1-2b9f-45fc-bfda-a525cdf8679f
# ╟─48db4cd0-15f0-11ee-0796-d3fffb0d18a2
# ╟─2a7ab774-12b0-4a75-8266-891d8411c1f2
# ╟─7dc5f83b-7a07-4741-9a2e-42e9246a8999
# ╟─c5d0dd04-241e-48a8-944d-dfa61edeface
