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

function G = grass_build_graph_hcp(ID,gtype,opts,overwrite)
% Constructs cerebral hemisphere cortex (CHC) brain graph using HCP data.
% The graph is specific to the left or right hemisphere. Voxels within mask
% become graph vertices. Graph edges defined base on two principles: (i)
% 26-neighborhood connectivity of voxels, (ii) pruning out anatomically
% invalid edges resulting from (i) using brain pial and white surface
% files. Pruning is time consuming since each individual edge derived based
% on (i) should be examined to infer whether or not it passes through the
% pial/white surface. This can be done more efficiently by paralleliing the
% procedure.
 
% Inputs: 
%   ID: HCP subject ID.
%   gtype: CHC graph type; gm*h.res****.spaceT1w
%   opts: structure with a bunch of necessary fields. 
% Output: 
%   G: graph structure with adjacency matrix, filenames, etc.  

%-Initiate G structure.
G = grass_get_G_hcp(ID,gtype,opts);

%-Check if G exists.
sts = checkexist(G,overwrite);
if isstruct(sts)
    G = sts;
    return;
end

% Build A.
A = buildA(G);

% Prune A.
G = pruneA(A,G,opts);

%-Save G.
save(G.fname,'G','-v7.3');
end

%==========================================================================
function sts = checkexist(G,overwrite)
chk = false;
if isstruct(overwrite)
    if ~overwrite.graph
        chk = true;
    end
elseif ~overwrite
    chk = true;
end
sts = [];
if chk
    if exist(G.f.G,'file')
        d = load(G.f.G);
        if isfield(d.G,'A')
            fprintf('\n..Graph already constructed.\n');
            sts = d.G;
            return;
        end
    end
end
end

%==========================================================================
function [A,indices,mask] = buildA(G)
% Based on:
% https://github.com/aitchbi/gwSPM/blob/master/gwspm_compute_adjacency.m

fprintf('\n..Determining voxel adjacencies.. ')

%-Extract cerebral hemisphere cortex mask.
mask = getmask(G);
    
dim = size(mask);

N = numel(mask); 

indices = find(mask);

[ii,jj,kk] = ind2sub(dim,indices);

%-Determine voxel adjacencies.
ci = repmat(ii,13,1);
cj = repmat(jj,13,1);
ck = repmat(kk,13,1);

ni = [ii;ii+1;ii+1;ii+1;ii;ii;ii;ii+1;ii+1;ii+1;ii+1;ii+1;ii+1];
nj = [jj+1;jj;jj-1;jj+1;jj;jj+1;jj+1;jj;jj;jj-1;jj+1;jj+1;jj-1];
nk = [kk;kk;kk;kk;kk+1;kk-1;kk+1;kk-1;kk+1;kk-1;kk-1;kk+1;kk+1];

maskZ = cat(2,mask,zeros(dim(1),1,dim(3)));
maskZ = cat(1,maskZ,zeros(1,dim(2)+1,dim(3)));
maskZ = cat(3,maskZ,zeros(dim(1)+1,dim(2)+1,1));

valid = (ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2) & nk>=1 & nk<=dim(3));
ni = ni(valid);
nj = nj(valid);
nk = nk(valid);
ci = ci(valid);
cj = cj(valid);
ck = ck(valid);

d = sub2ind(size(maskZ),ni,nj,nk);
ee = maskZ(d);
valid = logical(ee);
ni = ni(valid);
nj = nj(valid);
nk = nk(valid);
ci = ci(valid);
cj = cj(valid);
ck = ck(valid);

cInd = sub2ind(dim,ci,cj,ck);
nInd = sub2ind(dim,ni,nj,nk);

%-Build A.
A = sparse([cInd,nInd],[nInd,cInd],ones(1,2*numel(ni)),N,N);

%-Remove empty rows & cols.
c = find(~sum(A,1));
c = c(~ismember(c,indices));
A(:,c) = [];
A(c,:) = [];

fprintf('\n..done.\n')
end

%==========================================================================
function v = getmask(G)
% Extracting tissue mask, and then resample/reslice if necessary.
% The input file is the T1w/ribbon.nii reulting from FreeSurfer.

f_o =  G.f.mask; 

f_tmp = fullfile(fileparts(f_o),'temp.nii');

%-Extract tissue mask.
h_rb = spm_vol(hb_gunzip(G.f.source));
v = spm_read_vols(h_rb);
switch G.tissue
    case 'gmlh'
        t = 3;
    case 'gmrh'
        t = 42;
end

v = ismember(v,t);
h0 = struct();
h0.fname = f_tmp;
h0.dim = h_rb.dim;
h0.dt = [spm_type('uint8') h_rb.dt(2)]; % See NOTE 1
h0.mat = h_rb.mat;
h0.pinfo = [h_rb.pinfo(1:2);0]; % See NOTE 2
h0 = spm_create_vol(h0);
spm_write_vol(h0,v);
clear v;

%-Resample/Reslice.
interp = 1;
res_o = str2double(G.res)/1e3;
d = abs(diag(h_rb.mat));
if isequal(d(1),d(2),d(3))
    res_i = num2str(d(1));
else
    res_i = 'n/a';
end
if res_o>=res_i
    apprch = 'approach1';
else
    apprch = 'approach2';
end
if ~strcmp(res_i,num2str(res_o))
    f_tmp = hb_resample_vol(f_tmp,res_o,interp,f_tmp,apprch);
end

%-Threshold & Binarize.
% -If f_tmp is not resampled/resliced, it's already binary.
% -If f_tmp is resampled/resliced, it's aready thresholded & binarized,
%  since we set its header h.dt=2, i.e., 'uint8' type; See NOTE 3.

%-Make connected.
h = spm_vol(f_tmp);
[v,d] = hb_make_connected(h,26);
fprintf('\n  Mask made connected:');
fprintf('\n  number of voxels removed: %d',d); 
h.fname = f_o; %G.f.mask;
spm_write_vol(h,v);

delete(f_tmp);
end

%==========================================================================
function G = pruneA(A0,G,opts)
% Remove edges that pass through pial surfaces, and then, remove nodes 
% that become disconnected as a result of this pruning.

if opts.parallelize_prune
    gcp;
    prl = true;
else
    prl = false;
end

fprintf('\n..Pruning anatomically invalid graph edges.. ');

h_mask = spm_vol(G.f.mask);
mask = spm_read_vols(h_mask);
inds_i = find(mask);

%-Prune edges passing through pial.
d = ml_batch_gifti(G.f.surface.pial,h_mask.mat);
Ap_pial = ml_prune_adjacency(A0,mask,d,'parallelize',prl);
d1 = length(find(A0));
d2 = length(find(Ap_pial));
fprintf('\n  %d edges removed. [based on pial surface]',(d1-d2)/2);

%-Prune edges passing through white.
d = ml_batch_gifti(G.f.surface.white,h_mask.mat);
Ap_white = ml_prune_adjacency(Ap_pial,mask,d,'parallelize',prl);
d1 = length(find(Ap_pial));
d2 = length(find(Ap_white));
fprintf('\n  %d edges removed. [based on white surface]',(d1-d2)/2);

Ap = Ap_white;

%-Clean A.
[A,indA,mask] = ml_clean_adjacency(Ap,1,mask);
d1 = length(inds_i);
d2 = length(find(indA));
fprintf('\n  %d nodes removed.',d1-d2);
indices = inds_i(indA);
if setdiff(indices,find(mask)), error('fishy.'); end
chk = setdiff(inds_i,indices); %removed nodes
if ~isequal(chk,inds_i(~indA)),error('fishy.'); end

%-Build volume highlighting removed voxels from G.f.mask.
h = h_mask;
v = spm_read_vols(h);
v(inds_i(~indA)) = 2;
[p,n,e] = fileparts(h.fname); 
h.fname = fullfile(p,[n,'.postpruning_cleaned_nodes_marked',e]); 
spm_write_vol(h,v);
gzip(h.fname);
delete(h.fname);

%-Update G.
G.N = size(A,1);
G.A = A;
G.indices = indices;
G.pruning.A_diff_pre_post_pial_pruning = A0-Ap_pial; 
G.pruning.A_diff_pre_post_white_pruning = Ap_pial-Ap_white; 
G.pruning.A_removed_post_pruning = Ap(~indA,~indA);
G.pruning.ind_pre_pruning_A_remained_post_pruning = indA; 
G.pruning.indices_of_mask_cleaned_post_pruning = inds_i(~indA);
G.pruning.mask_with_post_pruning_cleaned_voxels_marked = h.fname;
G.pruning.info = fullfile(p,[n,'.pruning_info.txt']); 

%-Info file for plotting prunned edges.
fid = fopen(G.pruning.info,'wt');
fprintf(fid,'=========================================================\n');
fprintf(fid,'%%To plot prunned edges: \n');
fprintf(fid,'%%A0: pre-pruning A. \n');
fprintf(fid,'%%App: post-pial-pruning A. \n');
fprintf(fid,'%%Apw: post-white-pruning A. \n');
fprintf(fid,'\n');
fprintf(fid,'markCleanedVoxelsPostPruning = true; %%false \n');
fprintf(fid,'Adiffp = G.pruning.A_diff_pre_post_pial_pruning; \n');
fprintf(fid,'Adiffw = G.pruning.A_diff_pre_post_white_pruning; \n');
fprintf(fid,'Armd = G.pruning.A_removed_post_pruning; \n');
fprintf(fid,'indA = G.pruning.ind_pre_pruning_A_remained_post_pruning;\n');
fprintf(fid,'I = G.pruning.indices_of_mask_cleaned_post_pruning;\n');
fprintf(fid,'Apw = Adiffw-Adiffw;     %%empty matrix \n');
fprintf(fid,'Apw(indA,indA) = G.A;    %%place cleaned A \n');
fprintf(fid,'Apw(~indA,~indA) = Armd; %%put back removed rows&columns \n');
fprintf(fid,'App = Apw+Adiffw;        %%put back white-pruned edges \n');    
fprintf(fid,'A0 = App+Adiffp;         %%put back pial-pruned edges \n');  
fprintf(fid,'hb_gunzip(G.f.mask); \n');
fprintf(fid,'mask = spm_read_vols(spm_vol(G.f.mask)); \n');
fprintf(fid,'if markCleanedVoxelsPostPruning \n');
fprintf(fid,'    d = find(mask); \n');
fprintf(fid,'    mask(d) = 2; \n');
fprintf(fid,'    mask(I) = 1; \n');
fprintf(fid,'end \n');
fprintf(fid,'sl = 88;       %%axial slice number \n');
fprintf(fid,'d1 = ''grid''; %%''diff'',''none'' \n');
fprintf(fid,'d2 = true;     %%false \n');
fprintf(fid,'d3 = {''r'';''g''}; \n');
fprintf(fid,'figure; \n');
fprintf(fid,'A_pre = {A0;App}; \n');
fprintf(fid,'A_post = {App;Apw}; \n');
fprintf(fid,'hm_plot_adjacency_diff(A_pre,A_post,mask,sl,d1,d2,d3);\n');
fprintf(fid,'=========================================================\n');
fclose(fid);

% Compress.
gzip(G.f.mask);
delete(G.f.mask);

fprintf('\n..done.\n')
end

%==========================================================================
%-NOTES-
%
% NOTE 1 ------------------------------------------------------------------ 
% This is the data type --- See spm_type.m
% Here, we set it to spm_type('uint8'), i.e., 2, since our file is:
% -for gmlh/rh: a logocal file, storing only 2 numbers: 0 & 1.
% -for wmlh/rh: a file storing either 2 (0,1) or only 3 (0,1,2) numbers.
% These few values fall within the accepatable range for 'uint8': 0-2^8-1
% Using this small data type, significantly reduces the file size:
% 0.7 mm^3 res: 84 MB ('float32') > 0.7 mm^3 res: 21 MB ('unit8').
% For instance, the 1 mm^3 res file will become 7 MB ('uint8').
%
% NOTE 2 ------------------------------------------------------------------ 
% The third element of pinfo, refers to off set in bytes; for instance, see
% the following link for some further info:
% https://nipy.org/nibabel/devel/spm_use.html
% It is imprtant to change this, since we are saving each plane as a logical
% array and thus essentially there needs not be any offset in bytes.
% If we take the value from h_rb/mask.pinfo, then that can be problematic,
% since that value is potentially associated to the other data types that
% were used for saving those volumes, probably as a float32 data type.
% spm_read_vols, etc. will run into problem when wanting to load the files
% this property is not correctly set.
%
% NOTE 3 ------------------------------------------------------------------
% If the input volume is 'logical', like the volumes create here, then the
% output resampled volume is forced to be logical, so no need for
% thresholding. This is what spm_write_vol.m and spm_write_plane.m impose
% by seeing that the h.dt(1) is 2, i.e., 'uint8': unsigned integer, with
% min-max: 0-2^8-1; see NOTE 1. Thus, values in range [0,1] are rounded to
% closest integer value, i.e., 0 or 1.
