function [h,v] = ml_load_nifti(fname,volumes)
% ML_LOAD_NIFTI Extracts and loads a NIFTI-file.
%   ML_LOAD_NIFTI(fname)
%       Makes sure that the file fname exists. If not the function looks
%       for [fname '.gz'] and extracts it.
%   h = ML_LOAD_NIFTI(fname)
%       Loads the header of the NIFTI-file specified by fname.
%   [h,v] = ML_LOAD_NIFTI(fname)
%       Loads the header and volume data of the NIFTI-file specified by
%       fname.
%   ... = ML_LOAD_NIFTI(fname,volumes)
%       If the NIFTI-file is a 4D volume only the volumes specified in
%       volumes are loaded.
%
%   Author:
%       Martin Larsson
%       June 2017

% Improvements in efficiency 12/2017 D

%     ml_extract(fname);

    if nargout > 0
        
        if nargin > 1
            if ~iscolumn(volumes)
                volumes = volumes';
            end
            volumes = num2str(volumes);
            files = strcat(fname,',',volumes);
        else
            files = fname;
        end
        
        h = spm_vol(files);
%         h = cell2mat(h);
        
        if nargout > 1
            v = spm_read_vols(h);
        end
    end
end
