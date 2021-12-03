export QsosedModel
struct QsosedModel <: Model
    parameters::Parameters
    bh::BlackHole
    warm::Warm
    corona::Corona
end

function QsosedModel(parameters::Parameters)
    bh = BlackHole(parameters)
    corona = Corona(bh, parameters)
    warm = Warm(corona, parameters)
    return QsosedModel(parameters, bh, warm, corona)
end
function QsosedModel(config)
    parameters = Parameters(config)
    return QsosedModel(parameters)
end


