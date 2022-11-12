module CustomStructs
export SystemParams, LeadParams, LeadParamsSF

struct SystemParams
    N   :: Int64
    HS  :: Array
    En  :: Array
end    
    
mutable struct LeadParams
    muL :: Float64
    muR :: Float64
    gLu  :: Float64
    gLd  :: Float64
    gRu  :: Float64
    gRd  :: Float64
    tL  :: Float64
    tR  :: Float64
end

mutable struct LeadParamsSF
    muL :: Float64
    muR :: Float64
    gLu  :: Number
    gLd  :: Number
    gLud :: Number
    gLdu :: Number 
    gRu  :: Number
    gRd  :: Number
    gRud :: Number
    gRdu :: Number
    tL  :: Float64
    tR  :: Float64
end


end
