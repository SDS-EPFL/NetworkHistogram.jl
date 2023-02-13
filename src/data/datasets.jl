using HTTP, CodecZstd, TranscodingStreams

const url_ref = "https://networks.skewed.de"

include("utils.jl")

function get_netzschleuder_network(name::String)
    url = joinpath(url_ref, "net", name, "files", name * ".gt.zst")
    res = HTTP.get(url)

    if res.status != 200
        error("Error downloading network" * res.status)
    end

    decompressed = Base.IOBuffer(transcode(ZstdDecompressor, res.body))
    return readgt(decompressed)
end
