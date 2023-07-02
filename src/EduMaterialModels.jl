module EduMaterialModels
using Newton, Tensors, StaticArrays, ForwardDiff
using MaterialModelsBase
import MaterialModelsBase as MMB

include("utils/svoigt.jl")    # Conversion to static voigt/mandel #Tensors.jl#194
include("utils/utils.jl")     # Useful general functions 
include("utils/driver.jl")    # Simple driver to simulate a strain history

# Material behaviors
include("Elastic.jl")           # Elastic behavior, can also be used inside other materials
include("Plastic.jl")           # Rate independent and rate dependent plasticity
include("PlasticExplicit.jl")   # Special implementation of Forward Euler time integration
include("ViscoElastic.jl")      # 1D Zener material (uniaxial stress) 

end
