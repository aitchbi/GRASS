% Part of GRAph-based Spatial Smoothing (GRASS): 
% https://github.com/aitchbi/GRASS
%
% Hamid Behjat

function grass_get_gsigs(G,f_fmri,f_gsigs,overwrite)
% GRASS_GET_GSIGS does the following in (for HCP data):
% 1. transfers & resamples an fMRI file from MNI to T1w,
% 2. extracts graph signals.
%
% Inputs: 
%   G: graph structure with fields f, indices, etc. 
%   f_fmri: fMRI file; absolute address. 
%   f_gsigs; graph signals file to be written; absolute address.
%   overwrite: (optional, default: false) overwrite existing files? 

if ~exist('overwrite','var') || isempty(overwrite)
    overwrite_gsigs = false;
else
    if isstruct(overwrite)
        overwrite_gsigs = overwrite.gsigs;
    else
        overwrite_gsigs = overwrite;
    end
end

[p,n,e] = fileparts(f_gsigs);
f_gsigs_done = fullfile(p,[n,'_done',e]);

% Graph signals already extracted?
if exist(f_gsigs,'file')
    d1 = spm_vol(f_run);
    d2 = load(f_gsigs);
    chk = d1==size(d2.gsigs,2); % correct length?
    if chk==false
        delete(f_gsigs);
    end
else
    chk = false;
end

if all([chk,~overwrite_gsigs,exist(f_gsigs_done,'file')])
    fprintf('\n..Graph signals already extracted.')
    return;
end

%-Extract graph signals ---------------------------------------------------
h_fmri = spm_vol(f_fmri);
Nv = length(h_fmri);
inds = G.indices;
gsigs = zeros(length(inds),Nv);
[xx,yy,zz] = ind2sub(h_fmri(1).dim,inds);

for iV=1:Nv
    gsigs(:,iV) = spm_sample_vol(h_fmri(iV),xx,yy,zz,0);
end
save(f_gsigs,'gsigs','-v7.3');
done = true;
save(f_gsigs_done,'done');
end

%==========================================================================
function transf(source,dest)
if ~exist([source,'.gz'],'file')
    error('MNI fMRI volume missing.');
end
[~,~] = mkdir(fileparts(dest));
copyfile([source,'.gz'],[dest,'.gz'])
end

