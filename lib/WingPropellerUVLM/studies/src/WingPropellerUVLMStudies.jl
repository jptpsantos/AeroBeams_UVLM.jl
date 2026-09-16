module WingPropellerUVLMStudies
import WingPropellerUVLM
include("Convergence.jl")
include("MovingBlockDamping.jl")
include("moving_block_adapter.jl")
include("StudySettings.jl")
include("StudyRunner.jl")
include("ChangProblem.jl")
const LOADED_SOURCE=Convergence.SOURCE()

export AerodynamicConvergenceStudy, AeroelasticConvergenceStudy
export create_AerodynamicConvergenceStudy, create_AeroelasticConvergenceStudy
export study_settings, canonical_settings, legacy_settings, solve!, inspect_selection
export create_ChangModel, create_ChangDynamicProblem, ChangDynamicProblem
end

