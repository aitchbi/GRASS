% % [10 feb 2018] added option to input precomputed cheby coeffs
% 
% added version '2' option for showing the sum of abs values rather than 
% sum of squared of abs values. In particular, if the input 'g' was esigned
% as 'uniform_meyer_type_v2', this version should be used for ploting, to
% correctly obtain the partition of unity constraint. Note that I use this 
% _v2 design for the defining subdomains for DILi and benchmarking.   
% [25 June 2017]
%
% [21 June 2017]
% Martin added a fix so it can handle the color profile of an arbitrary 
% number of kernels.   
%
%
% [21 March 2017]
% added extra checks.  If G is a structure, it may still not have the
% field 'E', since we might have skipped computing it.
%
%


%SPG_VIEW_DESIGN Illustrate a system of spectral kernels.
%
% This function is an extended version of the function:
%
% sgwt_filter_design
%
% from the the SGWT toolbox.
%
% To see examples of usage and the required inputs, please refer to
% either of the following demo functions:
%
% spg_demo_ (minnesota | alameda |cerebellum)
% spg_demo_uniformMeyerType
% spg_demo_ (your_data_frame | your_data_decompose)
%
% ------------------------------------------------------------
% Copyright (C) 2016, Hamid Behjat.
% This file is part of SPG (signal processing on graphs) package - v1.00
%
% Download: miplab.epfl.ch/software/
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [f,x,maxY,maxY_cheby] = spgg_view_design(g,arange,varargin)

control_params={...
    'Graph',0,...
    'eigsToInterp',[],...
    'subSampleWarping',0,...
    'warping','none',...
    'plotLineWidth',2,...
    'lambda',[],...
    'guiHandle',[],...
    'chebyOrder',[],...
    'onlyChebyApprox','no',...
    'plotOnlyLP',0,...
    'vers',1,... %[25 June 2017]
    'c',[],... % [10 feb 2018]
    'kernelCenters',[],... % Sep 2019. 
    'tightFrameBounds',[],... % Sep 2019. 
    'onlyPlotChebyApproxTightFrameDeviation',false}; % Oct 2019 - skip plotting individual kernel cheby ord approximations; just plot the resulting tight frame bound deviation. 

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

cOrd = chebyOrder;
GG = Graph;

if any([isa(GG,'struct'),isa(GG,'matlab.io.MatFile')]) ...
        && all(ismember({'N','E'},fieldnames(GG)))
    
    if numel(GG.E)==GG.N
        theEigs = GG.E;
        diffOfEigs = theEigs(2:end)-theEigs(1:end-1);
        indNonEqual = [1; find(logical(diffOfEigs)) + 1];
        noSpectrum = 0;
    else
        noSpectrum = 1;
    end
else
    noSpectrum = 1;
end


if ~noSpectrum
    if strcmp(warping,'none')
        if isempty(eigsToInterp)
            f = GG.E(indNonEqual);
            x = f;
        else
            f = eigsToInterp;
            x = f;
        end
    else
        if isempty(eigsToInterp)
            f = GG.E(indNonEqual);
            x = warping(indNonEqual)*GG.lmax;
        else
            if subSampleWarping
                dummy1 = GG.E(indNonEqual);
                interp_x = dummy1(1:subSampleWarping:end-1);
                interp_x(end+1) = dummy1(end);
                
                dummy2 = warping(indNonEqual)*GG.lmax;
                interp_y = dummy2(1:subSampleWarping:end-1);
                interp_y(end+1) = dummy2(end);
                
                f = eigsToInterp(eigsToInterp <= interp_x(end));
                x = abs(gsp_mono_cubic_warp_fn(interp_x,interp_y,f));
            else
                f = eigsToInterp;
                interp_x = GG.E(indNonEqual);
                interp_y = warping(indNonEqual)*GG.lmax;
                x = abs(gsp_mono_cubic_warp_fn(interp_x,interp_y,f));
            end
        end
    end
    
elseif ~isempty(eigsToInterp) && strcmp(warping,'none')
    f = eigsToInterp;
    x = f;
else
    x = linspace(arange(1),arange(2),1e3);
    f = x;
end

J = numel(g)-1;
co = get(gca,'ColorOrder');
G = zeros(size(x));
G_cheby = zeros(size(x));

if ~isempty(cOrd) 
    c = cell(size(g));
    if length(g)>1 && length(cOrd)==1
        cOrd = cOrd*ones(size(g));
    end
    for k = 1:numel(g)
        cord = cOrd(k);
        c{k} = sgwt_cheby_coeff(g{k},cord,cord+1,arange);
    end
end

for n=0:J
    color = co(mod(n,size(co,1))+1,:);
    if isempty(guiHandle)
        if isempty(cOrd) && isempty(c)
            plot(f,g{1+n}(x),'Color',color,'LineWidth',plotLineWidth);
        elseif strcmp(onlyChebyApprox,'yes')
            plot(f,sgwt_cheby_eval(x,c{1+n},arange),'Color',color,...
                'LineWidth',plotLineWidth);
        else
            plot(f,g{1+n}(x),'Color',color,'LineWidth',plotLineWidth);
            if n==0
                hold on
            end
            if ~onlyPlotChebyApproxTightFrameDeviation
                plot(f,sgwt_cheby_eval(x,c{1+n},arange),'-.',...
                    'Color',color,'LineWidth',plotLineWidth);
            end
        end
    else
        if isempty(cOrd)
            plot(f,g{1+n}(x),'Color',color,'LineWidth',plotLineWidth,...
                'Parent',guiHandle);
        elseif strcmp(onlyChebyApprox,'yes')
            plot(f,sgwt_cheby_eval(x,c{1+n},arange),'Color',color,...
                'LineWidth',plotLineWidth,'Parent',guiHandle);
        else
            plot(f,g{1+n}(x),'Color',color,'LineWidth',plotLineWidth,...
                'Parent',guiHandle);
            if n==0
                hold(guiHandle,'on')
            end
            plot(f,sgwt_cheby_eval(x,c{1+n},arange),'-.','Color',color,...
                'LineWidth',plotLineWidth,'Parent',guiHandle);
        end
    end
    
    if n==0 && plotOnlyLP
        break
    end
    
    if n==0
        if isempty(guiHandle)
            hold on
        else
            hold(guiHandle,'on')
        end
    end
    
    if isempty(cOrd)
        if vers==1
            d = g{1+n}(x).^2;
        elseif vers==2
            d = abs(g{1+n}(x));
        end
        d(isnan(d)) = 0;
        G = G+d;
    else
        if vers==1
            d1 = g{1+n}(x).^2;
            d2 = sgwt_cheby_eval(x,c{1+n},arange).^2;
        elseif vers==2
            d1 = abs(g{1+n}(x));
            d2 = abs(sgwt_cheby_eval(x,c{1+n},arange));
        end
        d1(isnan(d1)) = 0;
        G = G+d1;
        d2(isnan(d2)) = 0;
        G_cheby = G_cheby+d2;
    end
end

if n==0 && plotOnlyLP
    
else
    if isempty(guiHandle)
        if isempty(cOrd)
            plot(f,G,'k:','LineWidth',1);
        elseif strcmp(onlyChebyApprox,'yes')
            plot(f,G_cheby,'k:','LineWidth',1);
        else
            plot(f,G,'k:','LineWidth',1);
            plot(f,G_cheby,'k-.','LineWidth',1);
        end
    else
        if isempty(cOrd)
            plot(f,G,'k:','LineWidth',1,'Parent',guiHandle);
        elseif strcmp(onlyChebyApprox,'yes')
            plot(f,G_cheby,'k:','LineWidth',1,'Parent',guiHandle);
        else
            plot(f,G,'k:','LineWidth',1,'Parent',guiHandle);
            plot(f,G_cheby,'k-.','LineWidth',1,'Parent',guiHandle);
        end
    end
end

leglabels{1}='h';
for j=1:J
    leglabels{1+j}=sprintf('g_{%d}',j);
end
leglabels{J+2}='G';

if ~isempty(lambda)
    y=-.5*ones(size(lambda,1),1);
    stem(lambda,y, 'k');
    leglabels{J+3}='Lambda';
end

maxY = max(G(:));
if isempty(cOrd)
    maxY_cheby = [];
else
    maxY_cheby = max(G_cheby(:));
end

if ~isempty(tightFrameBounds)
    d = 0:0.1:eigsToInterp(end);
    plot(d,(1+tightFrameBounds)*ones(length(d)),':r');
    plot(d,(1-tightFrameBounds)*ones(length(d)),':r');
end

if ~isempty(kernelCenters)
    d = 0:0.1:1.1;
    for n=1:length(kernelCenters)
        plot(kernelCenters(n)*ones(size(d)),d,':k','LineWidth',1)
    end
end

if isempty(guiHandle)
    hold off
else
    hold(guiHandle,'off')
end
end


function argselectAssign(variable_value_pairs)
% argselectAssign : Assign variables in calling workspace
%
%  function argselectAssign(variable_value_pairs)
%  
%  Inputs : 
%  variable_value_pairs is a cell list of form
%  'variable1',value1,'variable2',value2,...
%  This function assigns variable1=value1 ... etc in the *callers* workspace
%
%  This is used at beginning of function to simulate keyword argument
%  passing. Typical usage is
%
%  argselectAssign(control_params);
%  argselectCheck(control_params,varargin);
%  argselectAssign(varargin);
%
%  where control_params is a cell list of variable,value pairs containing
%  the default parameter values.
% 
% See also argselectCheck
%
% Author : David K. Hammond, EPFL LTS2
% Date : December, 2007
% Project : common utilities

% This file is part of the SGWT toolbox (Spectral Graph Wavelet Transform toolbox)
% Copyright (C) 2010, David K. Hammond. 
%
% The SGWT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWT toolbox.  If not, see <http://www.gnu.org/licenses/>.

for j =1:2:length(variable_value_pairs)
  pname=variable_value_pairs{j};
  pval=variable_value_pairs{j+1};
  assignin('caller',pname,pval);
end
end

function argselectCheck(control_params,varargin_in)
% argselectCheck : Check if control parameters are valid 
%
%  function argselectCheck(control_params,varargin_in)
%
%  Inputs: 
%  control_params and varargin_in are both cell arrays
%  that are lists of pairs 'name1',value1,'name2',value2,...
%
%  This function checks that every name in varargin_in is one of the names
%  in control_params. It will generate an error if not, otherwise it will 
%  do nothing
%
%  This is used at beginning of function to simulate keyword argument
%  passing. Typical usage is
%
%  argselectAssign(control_params);
%  argselectCheck(control_params,varargin);
%  argselectAssign(varargin);
%
%  where control_params is a cell list of variable,value pairs containing
%  the default parameter values.
%
% See also argselectAssign
%
% Author : David K. Hammond, EPFL LTS2
% Date : December, 2007
% Project : common utilities

% This file is part of the SGWT toolbox (Spectral Graph Wavelet Transform toolbox)
% Copyright (C) 2010, David K. Hammond. 
%
% The SGWT toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% The SGWT toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with the SGWT toolbox.  If not, see <http://www.gnu.org/licenses/>.

[param_names{1:(length(control_params)/2)}]=control_params{1:2:end};
% check that passed in control parameters are valid
for j=1:2:length(varargin_in)
  if(isempty(strmatch(varargin_in{j},param_names)))
    error(['Invalid control parameter : ',varargin_in{j},...
           ' Valid options are',sprintf('\n'),...
          sprintf('%s \n',param_names{:})]);
  end
end
end


