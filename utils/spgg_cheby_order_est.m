function [cOrds,e,G,Gc] = spgg_cheby_order_est(g,arange,opts)
% HB_CHEBY_ORDER_EST estimates necessary chebyshev polynomial order to
% estimate a kernel such that the maximum absolute value difference between
% the kernel and the chbyshev polynomial estimate is less than tol.
%
% INPUTS:
%   g: an Nx1 cell array of kernels.
%   arange: spectral range, typically [0,lmax].
%   opts.minOrds: an Nx1 vector of min cheby order to consider for each kernel.
%   opts.maxOrds: an Nx1 vector of max cheby order to consider for each kernel.
%   opts.tol:  
%       - tol.kernel: tolerence when checking absolute difference between
%       estimate an orginal kernel.
%       - tol.tightframe: tolerence when checking tight frame or partition
%       of unity condition. 
%   opts.sz: step size used in realizing the kernel and it's estimate.
%   opts.vers: 1 [default] or 2. Filter design version. 
%              1: sosks associated to tight frame.
%              2: sosks forming partition of unity. 
%              See spgg_filter_design.m for further details.
%
% OUTPUTS:
%    cOrds.kernel: an Nx1 vector specifying estimated necessary polynomial
%    order to estimat each kernel such that the error is less than that
%    specified by opts.tol.kernel.
%
%    cOrds.tightframe: an Nx1 vector specifying estimated necessary
%    polynomial order to estimat each kernel such that the tight frame
%    condiction is satisfied upto an error given by opts.tol.tighframe.
%
%    e.kernel: an Nx1 vector specifying max error for the estimated cheby estimate
%    of the kernels.
%
%    e.tightframe: cheby approximation error of sum of squares of kernels,
%    across spectrum.
%   yKernels: an Nx1 vector specifying estimated necessary order for each kernel.
%
%   G: tight frame condition across spectrum using exact kernels. 
%
%   Gc: tight frame condition across spectrum using approximated kernels. 
%
% Hamid Behjat
% Sep 2019 - Oct 2021.

if ~exist('opts','var') || isempty(opts)
    opts=struct;
end
if ~isfield(opts,'minOrds') || isempty(opts.minOrds)
    minOrds = ones(1,length(g));
else
    minOrds = opts.minOrds;
end
if ~isfield(opts,'maxOrds') || isempty(opts.maxOrds)
    maxOrds = 800*ones(1,length(g));
else
    maxOrds = opts.maxOrds;
end
if ~isfield(opts,'tol') || isempty(opts.tol)
    tol.kernel = 1e-4; % see NOTE 1. 
    tol.tightframe = 1e-4;
else
    tol = opts.tol;
end
if ~isfield(opts,'vers') || isempty(opts.vers)
    vers = 1; 
else
    vers = opts.vers;
    assert(any(vers==[1 2]));
end
if ~isfield(opts,'sz') || isempty(opts.sz)
    sz = 1e-4;
else
    sz = opts.sz;
end
if ~isfield(opts,'justCheckKernels') || isempty(opts.justCheckKernels)
    justCheckKernels = false;
else
    justCheckKernels = opts.justCheckKernels;
end
if ~isfield(opts,'parallelize') || isempty(opts.parallelize)
    parallelize = true;
else
    parallelize = opts.parallelize;
end

N = length(g);
x = arange(1):sz:arange(2);
y = zeros(size(g));
z = zeros(size(g));
tolk = tol.kernel;
if ~justCheckKernels
    tolf = tol.tightframe;
end

% Estimate cheby order for kernels ----------------------------------------
% [30-april-2021]
% parfor n = 1:N 
%     d1 = g{n}(x);
%     for o=minOrds(n):maxOrds(n)
%         d2 = pvt1(g{n},o,arange,x);
%         d3 = abs(d1(:)-d2(:));
%         if all(d3<tolk)
%             break;
%         end
%     end
%     d = 'tol.kernel too strict; either increase tol.kernel or maxOrd.';
%     assert(all(d3<tolk),d);
%     y(n) = o;
%     z(n) = max(d3);
% end

if parallelize 
    parfor n = 1:N
        [y(n),z(n)] = pvt3(g,n,x,minOrds,maxOrds,arange,tolk);
    end
else
    for n = 1:N
        [y(n),z(n)] = pvt3(g,n,x,minOrds,maxOrds,arange,tolk);
    end
end

e = struct;
e.kernel_initial = z;
cOrds = struct;
cOrds.kernel = y;

if justCheckKernels
    e = 'n/a';
    G = 'n/a';
    Gc = 'n/a';
    return;
end

% tight frame [vers=1] or partition of unity [vers=2] tolerance -----------
while true
    [d1,G,Gc] = pvt2(x,N,g,y,arange,vers);
    if all(d1<tol.tightframe)
        e.tightframe = d1;
        cOrds.tightframe = y;
        for n = 1:N
            d1 = g{n}(x);
            d2 = pvt1(g{n},y(n),arange,x);
            d3 = abs(d1(:)-d2(:));
            e.kernel_final = max(d3);
        end
        return;
    end
    I = find(d1>=tolf);
    I = I(randperm(length(I))); % enables faster search across kernels in for loop below. 
    d2 = false(size(y));
    nn = 1:N;
    for iI=1:length(I)
        d3 = x(I(iI));
        for n=nn(:)' % should be row vector
            if g{n}(d3)>1e-3
                d2(n) = true;
            end
        end
        nn = find(~d2); % skip kernels already flagged.
    end
    y = y+1*d2;
    if any(y(:)>maxOrds(:))
        d = 'tol.tightframe seemingly too strict;';
        error([d,' either increase tol.tightframe or maxOrd.']);
    end
end
end

%==========================================================================
function y = pvt1(gn,o,arange,x)
% cheby polynomial estimate 
d = sgwt_cheby_coeff(gn,o,o+1,arange);
y = sgwt_cheby_eval(x,d,arange);
end

%==========================================================================
function [z,G,Gc] = pvt2(x,N,g,y,arange,vers)
% frame tight enough?
G = zeros(size(x));
Gc = zeros(size(x));
for n = 1:N 
    d1 = g{n}(x);
    d2 = pvt1(g{n},y(n),arange,x);
    switch vers
        case 1
            % frame
            d1 = d1.^2;
            d2 = d2.^2;
        case 2
            % partition of unity
    end
    G  = G+d1;
    Gc = Gc+d2;
end
z = abs(G-Gc);
end

%==========================================================================
function [o,e] = pvt3(g,n,x,minOrds,maxOrds,arange,tolk)
d1 = g{n}(x);
for o=minOrds(n):maxOrds(n)
    d2 = pvt1(g{n},o,arange,x);
    d3 = abs(d1(:)-d2(:));
    if all(d3<tolk)
        break;
    end
end
d = 'tol.kernel too strict; either increase tol.kernel or maxOrd.';
assert(all(d3<tolk),d);
e = max(d3);
end

% NOTE 1 
% Be cautious if using a tolrenace larger than 1e-4. If tol is too large,
% atoms will show spurios profiles, and filtered signals will be corrupted.
% For instance, on CHC graph atoms, a tol of 1e-2 leads to a Cheby order
% estimate that results in heat-kernel atoms with corrupted profiles around
% the center of the atom, and the resulting filtered volumes do not look
% nice.
