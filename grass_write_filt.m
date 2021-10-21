% Part of GRAph-based Spatial Smoothing (GRASS): 
% https://github.com/aitchbi/GRASS
%
% Hamid Behjat
% October 2021

function grass_write_filt(G,f_run,tau,coeffs)
[p,n,e] = fileparts(f_run);
f_graphfilt = fullfile(p,[n,'_graph_filt_',num2str(tau),e]);
h_run = spm_vol(f_run);
for iV = 1:length(h_run)
    if iV==1
        h = struct;
        h.fname = f_graphfilt;
        h.dim = h_run(1).dim;
        h.mat = h_run(1).mat;
        h.dt = h_run(1).dt;
        if exist(h.fname,'file')
            delete(h.fname)
        end
        h = spm_create_vol(h);
    end
    h.n(1) = iV;
    vol = zeros(h.dim);
    vol(G.indices) = coeffs(:,iV);
    spm_write_vol(h,vol);
end
end