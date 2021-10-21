clc
clear 

%-User settings.
dir_spm12 = '/address/to/spm12';
dir_hcp   = '/address/to/HCP_database';  % see NOTE 1
dir_save  = 'RESULTS';                   % see NOTE 2
gtype     = 'gmlh.res2000.spaceT1w';     % see NOTE 3
task      = 'tfMRI_EMOTION_LR';          % see NOTE 4
tau       = 10; % heat kernel parameter; value >0
ID        = 100307; % HCP subject ID
overwrite = false; % see NOTE 5

%-Preparotary stuff.
[opts,dir_save] = prepare(dir_hcp,dir_save,dir_spm12);

%-Build CHC graph.
grass_build_graph_hcp(ID,gtype,opts,overwrite);

%-Run GRASS.
grass_run_hcp(ID,gtype,task,tau,opts,overwrite); % see NOTE 6


%==========================================================================
% NOTE 1 ------------------------------------------------------------------ 
% This is the HCP root directory contating subject directories. Individual
% subject directory, which are in the format xxxxxx, where each x is a
% non-zero digit. The assumption is that you have correctly extracted the
% HCP data so that the data follows the HCP data strcuture as described in:
% https://www.humanconnectome.org/storage/app/media/documentation/s1200/...
% HCP_S1200_Release_Reference_Manual.pdf
% (see p.94, section: 'Directory structure for preprocessed MR data')
%
% NOTE 2 ------------------------------------------------------------------ 
% Results are written in a new directory, as specified by dir_save; specify
% either a name to create a directory in the same directory where this
% script is, or specify full address to a directory eleswhere. If you want
% to write results in same directory as your HCP data, set: 
% dir_save = dir_hcp;
% 
% NOTE 3 ------------------------------------------------------------------
% gtype specifies the cerebral hemisphere cortex (CHC) graph type. It is in
% this format: gm*h.res****.spaceT1w. Three settings are determined in this
% naming: (i) left or right hemisphere, (ii) resolution of the design, and
% (iii) space. res1000 means 1 mm^3, res1250 means 1.25 mm^3, and so on.
% Space can be either T1w or MNI, where the former means the acpc space as
% defined in the T1w folder of HCP data, and MNI means MNINonLinear folder
% of HCP data. Example options include:
%   'gmlh.res2000.spaceT1w'
%   'gmlh.res1750.spaceT1w'
%   'gmlh.res1500.spaceT1w'
%   'gmlh.res1250.spaceT1w'
% All above options can be created for right hemisphere by replacing gmlh
% with gmrh. 
%
% NOTE 4 ------------------------------------------------------------------
% Task can be any of the following from the HCP database:
%   'tfMRI_EMOTION_LR'
%   'tfMRI_GAMBLING_LR'
%   'tfMRI_LANGUAGE_LR'
%   'tfMRI_MOTOR_LR'
%   'tfMRI_RELATIONAL_LR'
%   'tfMRI_SOCIAL_LR'
%   'tfMRI_WM_LR'
%   'rfMRI_REST1_LR'
%   'rfMRI_REST2_LR'
%   'tfMRI_EMOTION_RL'
%   'tfMRI_GAMBLING_RL'
%   'tfMRI_LANGUAGE_RL'
%   'tfMRI_MOTOR_RL'
%   'tfMRI_RELATIONAL_RL'
%   'tfMRI_SOCIAL_RL'
%   'tfMRI_WM_RL'
%   'rfMRI_REST1_RL'
%   'rfMRI_REST2_RL'
%
% NOTE 5 ------------------------------------------------------------------
% 'overwrite' can be simply a boolean, true/false; in this case any exiting
% file that has been already generated will be overwritten/non-overwritten,
% respectively.
% 'overwrite' can be a structure, and if so, it should include all the
% following fields:
% - overwrite.graph
% - overwrite.gsigs
% - overwrite.mniToAcpcResampledFmri
% If you are not sure about some of the fields, just set them to true to be
% sure to overwrite already existing, previoulsy generated files; but not
% that this is not an efficinet strategy, since for instance, if you want
% to perform GRASS using different values of tau, you for instance don't
% need to regenerate the graph each time.
% 
% NOTE 6 ------------------------------------------------------------------
% In the current implementation, to create graphs that ideally match an
% individual's cortex, the design is done in the original subject space,
% i.e. HCP T1w (acpc) folder; this is what the tag '*.spaceT1w' refers to
% in gtype. A such, HCP preprocessed fMRI data are mapped back to the acpc
% space, and smoothing is performed in this space. Smoothed volumes are
% saved in directory <subj-ID>/T1w/Results/<task>. In a near future
% release, an option will be added to enable graph design and spatial
% smoothing in MNI space, through setting gtype to '*.spaceMNI'.


%==========================================================================
function [opts,dir_save] = prepare(dir_hcp,dir_save,dir_spm12)
d = fileparts(mfilename('fullpath'));
addpath(fullfile(d,'utils'));
addpath(dir_spm12); 
addpath(genpath(fullfile(d,'utils','prune')));
addpath(genpath(fullfile(d,'utils','external','sgwt_toolbox')));
if ~contains(dir_save,filesep)
    dir_save = fullfile(d,dir_save);
end
if ~exist(dir_save,'dir')
    mkdir(dir_save);
end
opts = struct;
opts.hcp_root = dir_hcp;
opts.hcpsave_root = dir_save;
opts.parallelize_prune = true; 
opts.deleteResampledVolume = true; % fMRI in MNI resliced to acpc
opts.deleteExtractedGraphSignal = true;
end
