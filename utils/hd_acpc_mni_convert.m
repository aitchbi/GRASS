function hd_acpc_mni_convert(f_s,f_d,f_r,f_o,varargin)
%HD_ACPC_MNI_convert Converts a nifti volume or series from ACPC to MNI
%space and vice versa using displacement field
%   Apply non-linear conversion using displacement fields to convert a
%   nifti volume or volume series between MNI and ACPC spaces.
%
% Inputs:
%   f_s: full path to nifti file to be converted. 
%
%   f_d: full path of nifti file containing the displacement field used in 
%   the space conversion. If f_s in ACPC and f_r in MNI, f_d is the 
%   ACPC-to-MNI displacement. If f_s in MNI and f_r in acpc, f_d is the 
%   MNI-to-ACPC displacement.
%
%   f_r: full path to file used as reference for the dimensions
%   of the output. f_o: Full path to file in which to save the results.
%
%   f_o: full path of the file that will be created and saved.
%
%   varargin{1} [whichVols]: indices of volumes to convert from 4D f_s.
%   NOTE: use this option carefully, since the save name is not adapted.
%
% David Abramian, June 2017
% Hamid Behjat, Nov 2018 - Oct 2021.


%-Code notation guide:
%--------------------------------------------------------------------------
% *_s: source file.
% *_r: reference file.
% *_d: displacement file.
% *_o: output file, i.e., resliced file in new space.
%--------------------------------------------------------------------------

if nargin<5
    whichVols = [];
else
    whichVols = varargin{1};
end

if nargin<6
    sS = true;
else
    sS = varargin{2};
end

if isempty(whichVols)
   whichVols = 1:length(spm_vol(f_s)); 
end

Nv = length(whichVols);

p = fileparts(f_o);
[~,~] = mkdir(p);

%-Displacement file.
%--------------------------------------------------------------------------
[h_d,v_d] = ml_load_nifti(f_d);
mat_d = h_d.mat;

%-Reference file.
%--------------------------------------------------------------------------
h_r = ml_load_nifti(f_r,1);
mat_r = h_r.mat;
matinv_r = abs(inv(mat_r));
matinv_r(1:3,4) = 0;

dim_o = h_r.dim;

res_matinv_r = diag(matinv_r(1:3,1:3))';
for iV = 1:Nv
    if sS, hb_progress(iV,Nv,'tag','  '); end
    
    % source file
    d = strcat(f_s,',',num2str(whichVols(iV)));
    h_s = spm_vol(d);
    v_s = spm_read_vols(h_s);
    mat_s = h_s.mat;

    % remove NaN values; interpolation fails: erosion around NaNs
    v_s(isnan(v_s)) = 0;
    
    % voxel coordinates of output file
    [yy_o,xx_o,zz_o] = meshgrid(1:dim_o(2),1:dim_o(1),1:dim_o(3));

    %-Step1: ref space >> disp map
    %----------------------------------------------------------------------
    % positions in disp map
    A = affine3d((mat_d\mat_r)');
    [yy_d,xx_d,zz_d] = transformPointsForward(A,yy_o,xx_o,zz_o);

    % interpolate disp map to find values at positions
    disp_xx = interp3(v_d(:,:,:,1),yy_d,xx_d,zz_d);
    disp_yy = interp3(v_d(:,:,:,2),yy_d,xx_d,zz_d);
    disp_zz = interp3(v_d(:,:,:,3),yy_d,xx_d,zz_d);
    
    %-Step2: disp map >> ref space
    %----------------------------------------------------------------------
    % positions in ref space
    xx_r = xx_o + res_matinv_r(1)*disp_xx;
    yy_r = yy_o + res_matinv_r(2)*disp_yy;
    zz_r = zz_o + res_matinv_r(3)*disp_zz;
    
    %-Step3: ref space >> source space
    %----------------------------------------------------------------------
    % positions in source space 
    A = affine3d((mat_s\mat_r)');
    [yy_s,xx_s,zz_s] = transformPointsForward(A,yy_r,xx_r,zz_r);
    
    % interpolate source file to find values at positions
    v_o = interp3(v_s,yy_s,xx_s,zz_s);
    v_o(isnan(v_o)) = 0;
    
    %-Save output.
    %----------------------------------------------------------------------
    if iV==1%~exist('h_o','var')
        h_o = struct;
        h_o.fname = f_o;
        h_o.dim   = h_r.dim;
        h_o.mat   = h_r.mat;
        h_o.dt    = h_s.dt;
        h_o = spm_create_vol(h_o);
    end
    h_o.n(1) = iV;
    if ~isempty(h_s.private.timing)
        h_o.private.timing = h_s.private.timing;
    end
    
    spm_write_vol(h_o,v_o);
end
end
