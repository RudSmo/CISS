module CustomStructs
export SystemParams, LeadParams, LeadParamsSF

struct SystemParams
    N   :: Int64
    HS  :: Array
    En  :: Array
end    
    
mutable struct LeadParams
    muLu :: Float64
    muLd :: Float64
    muRu :: Float64
    muRd :: Float64 
    gLu  :: Float64
    gLd  :: Float64
    gRu  :: Float64
    gRd  :: Float64
    tLu  :: Float64
    tLd  :: Float64
    tRu  :: Float64
    tRd  :: Float64
end

mutable struct LeadParamsSF
    muLu :: Float64
    muLd :: Float64
    muRu :: Float64
    muRd :: Float64 
    gLu  :: Number
    gLd  :: Number
    gLud :: Number
    gLdu :: Number 
    gRu  :: Number
    gRd  :: Number
    gRud :: Number
    gRdu :: Number
    tLu  :: Float64
    tLd  :: Float64
    tRu  :: Float64
    tRd  :: Float64
end


end