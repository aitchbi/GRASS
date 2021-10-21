function [vol_o,n,f] = hb_make_connected(vol_i,conn)
%HB_MAKE_CONNECTED make input mask (bw image. 2D or 3D) single connected.
% The largest component will be kept and the remainder of voxels/components  
% will be deleted. 
% 
% Inputs:
%   vol_i: bw 2D image or 3D volume. Or spm_vol header.
%   conn: 4 or 8 for 2D. 6, 18 or 26 for 3D. See bwconncomp.m for details.
%
% Outputs:
%   vol_o: bw image/volume made single connected. 
%   n: number of voxels removed to make mask single connected.
%   f: fraction of of voxels removed to make mask single connected.
%
% Dependencies: SPM toolbox.
%
% Hamid Behjat 
% April/Nov.2020

if isstruct(vol_i)
    vol_i = spm_read_vols(vol_i);
end

d1 = bwconncomp(vol_i,conn);

[d2,d3] = max(cellfun(@length,d1.PixelIdxList));

vol_o = zeros(size(vol_i));

vol_o(d1.PixelIdxList{d3}) = 1;

n = nnz(vol_i)-d2; % number of voxels removed
f = n/nnz(vol_i);  % fraction of voxels removed
if f>0.05, warning('A large number of voxels were removed!'); end

