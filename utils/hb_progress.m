function hb_progress(n,N,varargin)
% HB_PROGRESS display progress.
%
%
%
%
% Hamid Behjat 
% March 2017

control_params = {
    'init',0,...
    'l',numel(num2str(N)),...
    'tag',''};

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

if n==1 || init
    if isempty(tag)
        fprintf(['\n','progress: ']);
    else
        fprintf(['\n',tag]);
    end
else
    fprintf(repmat('\b',1,2*l+1),n);
end
eval(['fprintf(''%-',num2str(l),'d/%-',num2str(l),'d'',n,N)'])



%     dummy = '\b';
%     for j=1:2*d
%         dummy = [dummy,'\b']; %#ok<AGROW>
%     end