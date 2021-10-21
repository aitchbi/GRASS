% Part of GRAph-based Spatial Smoothing (GRASS): 
% https://github.com/aitchbi/GRASS
%
% Hamid Behjat
% October 2021

function [coeff,chebyOrds] = grass_filt(G,f_gsigs,g,chebyOrds)
% Spectral graph filteirng of graph signals via polynomial approximation.
%
% Inputs:
%   G: graph structure.
%   f_gsigs: .mat storing graph signals in it's columns. 
%   g: cell array of spectral kernels.
%   chebyOrds: (opt) structure with fields 'kernel' and 'tightframe' each
%   of length(g), specifying the Chebyev order to use for each kernel. If
%   not inut, will be estimated. 
% 
% Outputs: 
%   coeff: 
%   chebyOrds: 
% Dependencies: 
% .sgwt_toolbox
% .spgg_cheby_order_est.m

if length(g)==1
    g{2} = @(x) 1-g{1}(x); % see NOTE 1.
    appended = true;
end

% Load graph signals.
d = load(f_gsigs);
S = d.gsigs;
assert(size(S,1)==G.N,'G and G sigs mismatch in size.');

% Compute Laplacian.
if ~isfield(G,'L') || isempty(G.L)
    G.L = sgwt_laplacian(G.A,'opt','normalized');
    G.lmax = sgwt_rough_lmax(G.L);
end
if ~isfield(G,'lmax') || isempty(G.lmax)
    G.lmax = sgwt_rough_lmax(G.L);
end
arange = [0,G.lmax];

% Determine Chebyshev order; see NOTE 2.
if ~exist('chebyOrds','var') || isempty(chebyOrds)
    opts = struct;
    opts.justCheckKernels = true;
    opts.tol.kernel = 1e-4;
    opts.parallelize = false;
    [chebyOrds,~] = spgg_cheby_order_est(g,arange,opts);
else
    assert(isstruct(chebyOrds),'chebyOrds should be a structure.');
    assert(isfield(chebyOrds,'kernel'));
end

c = cell(1,length(g));
for k = 1:length(g)
    d = chebyOrds.kernel(k);
    c{k} = sgwt_cheby_coeff(g{k},d,d+1,arange);
end

% Decompose signals.
coeff = sgwt_cheby_op(S,G.L,c,arange);

if appended
    coeff = coeff{1};
end
end

% NOTE 1
% g{2} is not used; sgwt_cheby_op needs at least 2 kernels as input.
%
% NOTE 2
% Filtering is performed through polynomial approximation of the kernels. A
% Chebyshev polynomial of suitable order is estimated to approximate each
% spectral kernel. Then, a graph signal is filtered using the given kernel
% by applying the polynomial to the L matrix, and then applying the result
% to the graph signal; thus, only matrix operations on L are used,
% obviating the need to diagonalized L to obtain the eigenvectors.
