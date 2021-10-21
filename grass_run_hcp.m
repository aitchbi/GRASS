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

function grass_run_hcp(ID,gtype,task,tau,opts,overwrite)

% Get G and filenames.
G = grass_get_G_hcp(ID,gtype,opts,true);
f_run = G.f.fmri_graph.(task);
f_gsigs = G.f.signals.(task);

% Extract graph signals.
grass_get_gsigs_hcp(G,task,overwrite);

% GRASS.
fprintf('\n..Running GRASS.. ');
g{1} = @(x) exp(-x*tau); % heat kernel
gsigsfilt = grass_filt(G,f_gsigs,g);

% Write volumes & cleanup.
grass_write_filt(G,f_run,tau,gsigsfilt)
fprintf('done.\n');

% Cleanup.
if opts.deleteResampledVolume
    delete(f_run);
end
if opts.deleteExtractedGraphSignal
    delete(f_gsigs);
end
end

