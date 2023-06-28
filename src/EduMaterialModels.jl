module EduMaterialModels
using Newton, Tensors, StaticArrays, ForwardDiff
using MaterialModelsBase
import MaterialModelsBase as MMB

include("svoigt.jl")    # Conversion to static voigt/mandel #Tensors.jl#194
include("utils.jl")     # Useful general functions 
include("driver.jl")    # Simple driver to simulate a strain history

include("Elasticity.jl")
include("Plasticity.jl")
include("PlasticityExplicit.jl")

end
