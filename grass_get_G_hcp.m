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

function  G = grass_get_G_hcp(ID,gtype,opts,loadG)

if ~exist('loadG','var') || isempty(loadG)
    loadG = false;
end

if ~ischar(ID)
    ID = num2str(ID);
end

G = struct;
G.type = gtype;
G.subject = ID;
G.fname = [];
d = strfind(gtype,'.res');
G.tissue = gtype(1:d-1);
G.res    = gtype(d+(4:7));
if contains(gtype,'T1w')
    G.space = 'T1w';
end
G.resTag   = ['.res',G.res];
G.spaceTag = ['.space',G.space];
G.rsTag    = [G.resTag,G.spaceTag];
G.trsTag   = [G.tissue,G.rsTag];
G.settingsTag = '';
G.neighb = 3;
G.f = [];

f = struct; %file paths structure.

% Directories -------------------------------------------------------------
f.hcp_root = opts.hcp_root;
f.hcpsave_root = opts.hcpsave_root;

f.T1w = fullfile(opts.hcp_root,ID,'T1w');
f.MNI = fullfile(opts.hcp_root,ID,'MNINonLinear');
f.MNI_results      = fullfile(f.MNI,'Results');

f.T1w_save         = fullfile(opts.hcpsave_root,ID,'T1w');
f.MNI_save         = fullfile(opts.hcpsave_root,ID,'MNINonLinear');
f.T1w_results_save = fullfile(f.T1w_save,'Results');
f.MNI_results_save = fullfile(f.MNI_save,'Results');

f.graphmain        = fullfile(f.T1w_save,'graph'); 
f.graph            = fullfile(f.graphmain,strrep(gtype,'.','_'));

[~,~] = mkdir(f.T1w_save); % outputs set to prevent warnings
[~,~] = mkdir(f.MNI_save);
[~,~] = mkdir(f.MNI_results_save);
[~,~] = mkdir(f.graph);

% G.mat -------------------------------------------------------------------
f.G = fullfile(f.graph,['G.',G.type,'.mat']);

if loadG
    d = load(f.G);
    G = d.G;
    return;
end

% Volumetric files --------------------------------------------------------
% f.source is used for mask extraction.
% f.mask is the file used for G design, which has been:
% 1. Extracted from f.source.
% 2. Resampled/resliced to desired resolution/space.
% 3. Made connected.

n = [G.tissue,G.resTag];
f.source = fullfile(f.T1w,'ribbon.nii');
f.mask = fullfile(f.graphmain,[n,'.spaceT1w.nii']);

% Surface files -----------------------------------------------------------
d_surf = fullfile(opts.hcpsave_root,ID,'T1w',ID,'surf'); %see NOTE 1.
[~,~] = mkdir(d_surf);

f.surface.pial = {fullfile(d_surf,[gtype(3:4),'.pial.surf.gii'])};
f.surface.white = {fullfile(d_surf,[gtype(3:4),'.white.surf.gii'])};

fn = fieldnames(f.surface);
for iFN=1:length(fn)
    d = f.surface.(fn{iFN});
    for iF=1:length(d)
        if ~exist(d{iF},'file')
            d1 = fullfile(opts.hcp_root,ID,'T1w',ID,'surf');
            [~,d2,d3] = fileparts(d{iF});
            sts = copyfile(fullfile(d1,[d2,d3]),d{iF});
            if ~sts, error('[HB] problem copying surface file.'); end
        end
    end
end

% Mapping files -----------------------------------------------------------
f.xfms = fullfile(f.MNI,'xfms');
f.xfms_save = fullfile(f.MNI_save,'xfms');
[~,~] = mkdir(f.xfms_save);

% Displacement field for mapping from MNI -> ACPC
f.disp_mni2acpc = fullfile(f.xfms,'standard2acpc_dc.nii');

% Preprocessed fMRI volumes in MNI space ----------------------------------
% address to fMRI data that have been coregistered with graph.

taskSets{1} = [
    {'tfMRI_EMOTION_LR'   }
    {'tfMRI_GAMBLING_LR'  }
    {'tfMRI_LANGUAGE_LR'  }
    {'tfMRI_MOTOR_LR'     }
    {'tfMRI_RELATIONAL_LR'}
    {'tfMRI_SOCIAL_LR'    }
    {'tfMRI_WM_LR'        }
    {'rfMRI_REST1_LR'     }
    {'rfMRI_REST2_LR'     }];

taskSets{2} = [
    {'tfMRI_EMOTION_RL'   }
    {'tfMRI_GAMBLING_RL'  }
    {'tfMRI_LANGUAGE_RL'  }
    {'tfMRI_MOTOR_RL'     }
    {'tfMRI_RELATIONAL_RL'}
    {'tfMRI_SOCIAL_RL'    }
    {'tfMRI_WM_RL'        }
    {'rfMRI_REST1_RL'     }
    {'rfMRI_REST2_RL'     }];

for iS=1:2
    tasks = taskSets{iS}; 
    for t = 1:length(tasks)
        task = tasks{t};
        n = [task '.nii'];
        f.fmri_mni.(task)      = fullfile(f.MNI_results,task,n);
        f.fmri_mni_save.(task) = fullfile(f.MNI_results_save,task,n);
    end
end

% Resliced fMRI volumes to have 1-to-1 voxel macth with G.f.mask ----------
% fMRI data that have been coregistered with graph mask.
for iS=1:2
    tasks = taskSets{iS};
    for t = 1:length(tasks)
        task = tasks{t};
        n = [task,G.rsTag,'.nii'];
        f.fmri_graph.(task) = fullfile(f.T1w_results_save,task,n);
    end
end

% fMRI graph signals ------------------------------------------------------
% signals extracted from
% graph coregistered fMRI
% volumes.
for iS=1:2
    tasks = taskSets{iS};
    for t = 1:length(tasks)
        task = tasks{t};
        n = ['G.',G.type,'.signals_',task,'.mat'];
        f.signals.(task) = fullfile(f.graph,n);
    end
end

G.fname = fullfile(f.graph,['G.',G.type,'.mat']);
G.f = f;
end