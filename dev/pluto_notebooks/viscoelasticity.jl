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
function plot_response(;t_ramp=1.0, t_hold=10.0, E1=10e3, E2=20e3, η=1e3, ϵ=0.001, num_steps=100)
	tf(x) = SymmetricTensor{2,1}(tuple(x))
	ϵ_a = [tf(x) for x in range(0,ϵ,num_steps)]
	ϵ_b = [tf(ϵ) for _ in 2:num_steps]
	ϵ_vec = vcat(ϵ_a, ϵ_b)
	t = vcat(collect(range(0, t_ramp, num_steps)), range(t_ramp, t_ramp+t_hold, num_steps)[2:end])
	
	m = EduMaterialModels.Zener1D(;E1, E2, η)
	stress_state = UniaxialStress()
	σ = EduMaterialModels.simulate_response(m, stress_state, ϵ_vec, t)
	fig = CM.Figure()
	ax = CM.Axis(fig[1,1]; xlabel="ϵ₁₁ [%]", ylabel="σ₁₁ [MPa]")
	CM.lines!(ax, 100*first.(ϵ_vec), first.(σ))
	CM.xlims!(ax, 0, 100*ϵ)
	CM.ylims!(ax, 0, 1000)

	ax = CM.Axis(fig[1,1]; xlabel="t [s]", ylabel="σ₁₁ [MPa]")
	CM.lines!(ax, t, first.(σ))
	CM.xlims!(ax, 0, maximum(t))
	CM.ylims!(ax, 0, 1000)
	
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
	ts_slider = @bind tstar Slider([2.0^n for n in -2:6]; show_value=true, default=1.0)
	n_slider = @bind n Slider(1.0:0.25:4.0; show_value=true, default=1.0)
	t_slider = @bind cycle_time Slider([2.0^n for n in -1:10], show_value=true, default=4.0)
	N_slider = @bind num_steps Slider([2^n for n in 3:10]; show_value=true, default=64)
	
	md"""
	# Visco-plasticity modeling
	## Intended learning outcomes
	1. Understand how a visco-plastic model differs from a rate-independent plasticity model

	## Model theory
	The main difference between a rate-independent plasticity model and a visco-plastic material model, is that we no longer have the KKT-conditions. In particular, we allow stress-states outside the yield surface, i.e. ``\varPhi>0``. The plastic multiplier is then typically given on a form like

	``
	\dot{\lambda} = \frac{1}{t^*} \eta(\varPhi), \quad \eta(x) \geq 0
	``

	The function ``\eta`` is called the over-stress function, and there exist different versions in the literature. In this document, we will use one of versions that are denoted Norton overstress functions,

	``\eta(\varPhi) = \left[ \frac{\langle \varPhi \rangle}{Y_0}\right]^n, \quad \langle x \rangle = \mathrm{max}(0, x)``

	Where Y_0 is the initial yield limit. The viscous parts then introduce two new parameters, the characteristic time, ``t^*``, and the exponent, ``n``. Otherwise, the parameters and equations are the same as for the rate-independent material model. 
		
	## Simulation parameters
	This notebook simulates the uniaxial stress response for ``\epsilon_{11}=\pm 1 \%`` during one cycle. You can adjust the material parameters, ``E``, ``Y_0``, ``H_\mathrm{iso}``, ``\kappa_\infty``, ``H_\mathrm{kin}``, ``\beta_\infty``, ``t^*``, and ``n``. In addition, you can choose how many steps are taken for each quarter cycle, as well as the total time for one cycle. 
	
	| param | value        | unit |
	| ----- | ------------ | ---- |
	| E     | $(E_slider)  | GPa  |
	| Y0    | $(Y0_slider) | MPa  |
	| Hiso  | $(Hi_slider) | GPa  |
	| κ∞    | $(κ∞_slider) | MPa  |
	| Hkin  | $(Hk_slider) | GPa  |
	| β∞    | $(β∞_slider) | MPa  |
	| t*    | $(ts_slider) | s    |
	| n     | $(n_slider)  | -    |
	| cycle time  | $(t_slider)  | s    |
	| steps | $(N_slider)  | -    |  
	
	"""
end

# ╔═╡ c5d0dd04-241e-48a8-944d-dfa61edeface
plot_response()

# ╔═╡ Cell order:
# ╟─4b6d68c1-2b9f-45fc-bfda-a525cdf8679f
# ╟─48db4cd0-15f0-11ee-0796-d3fffb0d18a2
# ╠═2a7ab774-12b0-4a75-8266-891d8411c1f2
# ╟─7dc5f83b-7a07-4741-9a2e-42e9246a8999
# ╟─c5d0dd04-241e-48a8-944d-dfa61edeface
