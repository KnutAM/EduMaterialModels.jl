function simulate_response(m::AbstractMaterial, stress_state, strain_history, timevector)
    length(strain_history) == length(timevector)
    state = initial_material_state(m)
    
    σ, dσdϵ, state = material_response(stress_state, m, strain_history[2], state, timevector[2]-timevector[1])
    stress_history = zeros(typeof(σ), length(strain_history))
    stress_history[2] = σ
    for i in 3:length(strain_history)
        Δt = timevector[i]-timevector[i-1]
        σ, dσdϵ, state = material_response(stress_state, m, strain_history[i], state, Δt)
        stress_history[i] = σ
    end
    return stress_history
end

