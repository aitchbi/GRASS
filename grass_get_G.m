function G = grass_get_G(opts)
G = struct;
G.ID      = opts.ID;
G.type    = opts.gtype;
G.subject = opts.ID;
G.tissue  = opts.gtype(1:strfind(G.type,'.res')-1);
G.res     = opts.gtype(strfind(G.type,'.res')+(4:7));
G.resTag  = ['.res',G.res];
G.fname   = []; 
G.neighb  = 3;
G.f               = struct;
G.f.source        = opts.f_source;
G.f.surface.pial  = opts.f_surface.pial;
G.f.surface.white = opts.f_surface.white;
G.f.graphmain     = fullfile(opts.dir_save,G.ID,'graph');
G.f.graph         = fullfile(G.f.graphmain,strrep(G.type,'.','_'));
G.f.mask          = fullfile(G.f.graph,[G.tissue,G.resTag,'.nii']);
G.f.G             = fullfile(G.f.graph,['G.',G.type,'.mat']); % = G.fname
if not(exist(G.f.graphmain,'dir'))
    mkdir(G.f.graphmain);
end
if not(exist(G.f.graph,'dir'))
    mkdir(G.f.graph);
end
end