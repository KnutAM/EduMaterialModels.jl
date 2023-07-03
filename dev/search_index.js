var documenterSearchIndex = {"docs":
[{"location":"pluto_notebooks/plasticity_md/","page":"Plasticity","title":"Plasticity","text":"Click \"Edit or run this notebook\" to see interactive options.","category":"page"},{"location":"pluto_notebooks/plasticity_md/","page":"Plasticity","title":"Plasticity","text":"<iframe src=\"plasticity.html\" style=\"width: 100%; height: 2500px; border: none> </iframe> ","category":"page"},{"location":"models/","page":"Material Models","title":"Material Models","text":"CurrentModule = EduMaterialModels","category":"page"},{"location":"models/#Material-models","page":"Material Models","title":"Material models","text":"","category":"section"},{"location":"models/","page":"Material Models","title":"Material Models","text":"EduMaterialModels.LinearIsotropicElasticity\nEduMaterialModels.J2Plasticity\nEduMaterialModels.J2ViscoPlasticity\nEduMaterialModels.Zener1D","category":"page"},{"location":"models/#EduMaterialModels.LinearIsotropicElasticity","page":"Material Models","title":"EduMaterialModels.LinearIsotropicElasticity","text":"LinearIsotropicElasticity(;kwargs...)\n\nConstruct a linear, isotropic, elastic material.  The following pairs of properties can be given as keyword arguments:\n\nG, K: Shear and bulk modulus \nE, G: Young's and shear modulus \nE, ν: Young's modulus and Poisson's ratio \nE, K: Young's and bulk modulus\n\n\n\n\n\n","category":"type"},{"location":"models/#EduMaterialModels.J2Plasticity","page":"Material Models","title":"EduMaterialModels.J2Plasticity","text":"J2Plasticity(;e, Y0, Hiso, κ∞, Hkin, β∞)\n\nThe J2Plasticity material model assumes linear elasticity combined  with rate independent von mises plasticity. Nonlinear isotropic (Voce) and kinematic (Armstrong-Frederick) hardening  are included as well.\n\n\n\n\n\n","category":"type"},{"location":"models/#EduMaterialModels.J2ViscoPlasticity","page":"Material Models","title":"EduMaterialModels.J2ViscoPlasticity","text":"J2ViscoPlasticity(;e, Y0, Hiso, κ∞, Hkin, β∞, n, tstar)\n\nThe J2ViscoPlasticity material model assumes linear elasticity combined  with rate dependent (Norton-type) von mises plasticity. Nonlinear isotropic (Voce) and kinematic (Armstrong-Frederick) hardening  are included as well. Specifically, the overstress function is \n\neta(Phi) = left fraclangle Phi rangleY_0right^n\n\n\n\n\n\n","category":"type"},{"location":"models/#EduMaterialModels.Zener1D","page":"Material Models","title":"EduMaterialModels.Zener1D","text":"Zener1D(;E1, E2, η)\n\nSpecial implementation of the Zener material model, which only works  for uniaxial stress, i.e. when stress_state=MaterialModelsBase.UniaxialStress() is given to material_response. \n\n\n\n\n\n","category":"type"},{"location":"pluto_notebooks/viscoplasticity_md/","page":"Viscoplasticity","title":"Viscoplasticity","text":"Click \"Edit or run this notebook\" to see interactive options.","category":"page"},{"location":"pluto_notebooks/viscoplasticity_md/","page":"Viscoplasticity","title":"Viscoplasticity","text":"<iframe src=\"viscoplasticity.html\" style=\"width: 100%; height: 2000px; border: none> </iframe> ","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = EduMaterialModels","category":"page"},{"location":"#EduMaterialModels","page":"Home","title":"EduMaterialModels","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package contains basic material models for educational purposes following the MaterialModelsBase interface. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Short doc-strings are available for each material model,  and other utility functions for running simulations.  However, the main purpose is to demonstrate the different material behaviors using  Pluto notebooks:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Plasticity\nViscoplasticity\nViscoelasticity","category":"page"},{"location":"utils/","page":"Utility functions","title":"Utility functions","text":"CurrentModule = EduMaterialModels","category":"page"},{"location":"utils/#Utility-functions","page":"Utility functions","title":"Utility functions","text":"","category":"section"},{"location":"utils/","page":"Utility functions","title":"Utility functions","text":"EduMaterialModels.simulate_response","category":"page"},{"location":"utils/#EduMaterialModels.simulate_response","page":"Utility functions","title":"EduMaterialModels.simulate_response","text":"A simple driver to simulate a response for a given material, m, a stress_state<:MaterialModelsBase.AbstractStressState, a strain history  strain_history::Vector{<:AbstractTensor}, and time history,  timevector::Vector{<:Number}. \n\n\n\n\n\n","category":"function"},{"location":"pluto_notebooks/viscoelasticity_md/","page":"Viscoelasticity","title":"Viscoelasticity","text":"Click \"Edit or run this notebook\" to see interactive options.","category":"page"},{"location":"pluto_notebooks/viscoelasticity_md/","page":"Viscoelasticity","title":"Viscoelasticity","text":"<iframe src=\"viscoelasticity.html\" style=\"width: 100%; height: 2500px; border: none> </iframe> ","category":"page"}]
}