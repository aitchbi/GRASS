% Part of GRAph-based Spatial Smoothing (GRASS): 
% https://github.com/aitchbi/GRASS
%
% This implementation of GRASS is tailored to the Human Connectome Project
% (HCP) database. The assupmtion is that the HCP data are extracted in such
% way that the file directories follows that described in:
% https://www.humanconnectome.org/storage/app/media/documentation/s1200...
% /HCP_S1200_Release_Reference_Manual.pdf; see p.94, section: 'Directory
% structure for preprocessed MR data'.
%
% Hamid Behjat
% October 2021

function grass_get_gsigs_hcp(G,task,overwrite,verbose)
% GRASS_GET_GSIGS does the following in (for HCP data):
% 1. transfers & resamples an fMRI file from MNI to T1w,
% 2. extracts graph signals.
%
% Inputs: 
%   G: graph structure with fields f, indices, etc. 
%   task: HCP task name, e.g. 'tfMRI_EMOTION_LR'. 
%   overwrite: (optional, default: false) overwrite existing files? 
%   opts: (optional) structure with fields runPar,hcpsave_root, etc. 

if ~exist('overwrite','var') || isempty(overwrite)
    overwrite_gsigs = false;
    overwrite_reslicedFmri = false;
else
    if isstruct(overwrite)
        overwrite_gsigs = overwrite.gsigs;
        overwrite_reslicedFmri = overwrite.mniToAcpcResampledFmri;
    else
        overwrite_gsigs = overwrite;
        overwrite_reslicedFmri = overwrite;
    end
end
if ~exist('verbose','var') || isempty(verbose)
    verbose = true;
end

f_gsigs = G.f.signals.(task);
[p,n,e] = fileparts(f_gsigs);
f_gsigs_done = fullfile(p,[n,'_done',e]);

% Graph signals already extracted?
if exist(f_gsigs,'file')
    d1 = hb_get_hcp_task_length(task);
    d2 = load(f_gsigs);
    chk = d1==size(d2.gsigs,2); % extracted siganls of length = 4D volume?
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

%-Reslice fMRI ------------------------------------------------------------
f_fmri_mni      = G.f.fmri_mni.(task);   
f_fmri_mni_save = G.f.fmri_mni_save.(task);
f_fmri_graph    = G.f.fmri_graph.(task); 

if exist(f_fmri_graph,'file') && overwrite_gsigs
    delete(f_fmri_graph);
end

if exist(f_fmri_graph,'file')
    d1 = hb_get_hcp_task_length(task);
    d2 = length(spm_vol(f_fmri_graph));
    chk = d1==d2;
    if chk==false
        delete(f_fmri_graph);
    end
else
    chk = false;
end

if ~overwrite_reslicedFmri && chk
    fprintf('\n..Resliced fMRI volumes exist.');
else
    f_d = hb_gunzip(G.f.disp_mni2acpc,G.f.xfms_save);
    
    if all(...
            ~[exist(f_fmri_mni_save,'file'),...
            exist([f_fmri_mni_save,'.gz'],'file')])
        transf(f_fmri_mni,f_fmri_mni_save);
    end
    
    hb_gunzip(f_fmri_mni_save);
    
    hb_gunzip(G.f.mask);
    
    % reslice; mni >> G.f.tissue
    fprintf('\n..Reslicing fMRI volumes (MNI to ACPC)..')
    fprintf('\n   MNI res: 2 mm cubic');
    d = str2double(G.res)/1000;
    if rem(d,1)
        fprintf('\n  ACPC res: %.02f mm cubic',d);
    else
        fprintf('\n  ACPC res: %d mm cubic',d);
    end
    hd_acpc_mni_convert(...
        f_fmri_mni_save,...
        f_d,...
        G.f.mask,...
        f_fmri_graph,...
        [],...
        verbose);
    fprintf('\n..done.\n');
    delete(f_fmri_mni_save); % .gz version remains
    delete(G.f.mask); % .gz version remains
    if ~isequal(G.f.hcp_root,G.f.hcpsave_root)
        delete([f_fmri_mni_save,'.gz']); % original .gz remains in hcp_root
    end
end

%-Extract graph signals ---------------------------------------------------
h_fmri = spm_vol(f_fmri_graph);
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

