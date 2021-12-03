export QsosedModel
struct QsosedModel <: Model
    parameters::Parameters
    bh::BlackHole
    warm::Warm
    corona::Corona
end

function QsosedModel(config_file::String)
    parameters = Parameters(config_file)
    bh = BlackHole(parameters)
    corona = Corona(bh, parameters)
    warm = Warm(corona, parameters)
    return QsosedModel(parameters, bh, warm, corona)
end


