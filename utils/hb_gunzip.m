function f_o = hb_gunzip(f_i,varargin)
% HM_GUNZIP attempts to unzip [f_i,'.gz'] if f_i does not exist. 
%
% Inputs:
%   f_i: nifti file; full path.
%   varargin{1} [writeroot]: where to extract [f_i,'.gz'] to. 
%   varargin{2} [overwrite]: overwrite existing unziped file?
%
% Outputs:
%   f_o: full path to unziped f_i.
%
% Hamid Behjat
% April/Sep 2020

if nargin==3
    overwrite = varargin{2};
else
    overwrite = false;
end

[~,~,e] = fileparts(f_i);
if strcmp(e,'.gz')
    f_i = f_i(1:end-length(e));
end

if nargin<2 || isempty(varargin{1})
    writeroot = [];
    f_o = f_i;
else
    writeroot = varargin{1};
    [d1,d2,d3] = fileparts(f_i);
    if strcmp(d1,writeroot)
        writeroot = [];
        f_o = f_i;
    else
        f_o = fullfile(writeroot,[d2,d3]);
    end
end

f_igz = [f_i,'.gz'];
if ~exist(f_o,'file') || overwrite 
    if ~exist(f_igz,'file')
        error(['File does not exist: ',f_igz]);
    end
    if isempty(writeroot)
        d = gunzip(f_igz);
    else
        d = gunzip(f_igz,writeroot);
    end
    assert(strcmp(f_o,d{1}),'[HB] fishy.');
end
end
