function GammaReduced = partial_trace(Gamma, modes_to_keep);



Index = sort([2*modes_to_keep-1 2*modes_to_keep]);
GammaReduced = Gamma(Index, Index);


