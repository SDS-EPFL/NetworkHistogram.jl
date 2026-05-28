module MakieExt

using NetworkHistogram
using Makie

Makie.convert_single_argument(A::SymArray) = Matrix(A)

end
