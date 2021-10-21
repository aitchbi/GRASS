function [Ac,Aind,maskc] = ml_clean_adjacency(A,comps,mask)
% ML_CLEAN_ADJACENCY Limit number of connected components.
%   [Ac,Aind] = ML_CLEAN_ADJACENCY(A,comps)
%       Limits the number of connected components in the adjacency matrix A
%       by only keeping the comps largest ones. Ac is the cleaned adjacency
%       matrix and Aind is a logical vector indicating what nodes in A is
%       kept in Ac, i.e. Ac=A(Aind,Aind).
%   [Ac,Aind,maskc] = ML_CLEAN_ADJACENCY(A,comps,mask)
%       mask is a tensor whose non-zero elements, indexed linearly,
%       correspond to the nodes in A. maskc is a copy of mask but with the
%       removed nodes set to zero, thus Ac and maskc share the same
%       relationship as A and mask do.
%
%   Author:
%       Martin Larsson
%       March 2017

    bins = conncomp(graph(A));
    bin_count = max(bins);
    
    % There is no cleaning to be done.
    if bin_count <= comps
        Ac = A;
        Aind = true(1,size(A,1));
        if nargout >= 3
            maskc = mask;
        end
        return;
    end

    % Sort connected componens based on size.
    bins_size = histcounts(bins,1:bin_count+1);
    [~,bin_ind] = sort(bins_size,'descend');

    % Keep the comps largest connected components.
    Aind = ismember(bins,bin_ind(1:comps));
    Ac = A(Aind,Aind);

    % Clean mask.
    if nargout >= 3
        indices = find(mask);
        maskc = mask;
        maskc(indices(~Aind)) = 0;
    end
end

