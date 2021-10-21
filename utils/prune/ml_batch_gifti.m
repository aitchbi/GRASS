function surf = ml_batch_gifti(file_names,mat)
% ML_BATCH_GIFTI Load GIFTI files as single mesh.
%   surf = ML_BATCH_GIFTI(file_names)
%       Loads multiple GIFTI files into a single mesh. The vertices will be
%       transformed into ACPC space.
%   surf = ML_BATCH_GIFTI(file_names,mat)
%       mat will be used to transform the vertices from ACPC space such
%       that newV=mat\V.
%
%   Author:
%       Martin Larsson
%       March 2017

    if nargin < 2
        mat = eye(4);
    end

    n = length(file_names);
    
    surfs = cell(n,1);
    face_count = 0;
    vert_count = 0;
    
    for i=1:n
        surfs{i} = gifti(file_names{i});
        face_count = face_count + size(surfs{i}.faces,1);
        vert_count = vert_count + size(surfs{i}.vertices,1);
        surfs{i}.vertices = ml_transform_vertices(mat\surfs{i}.mat,...
            surfs{i}.vertices);
    end
    
    surf.faces = zeros(face_count,3,class(surfs{i}.faces));
    surf.vertices = zeros(vert_count,3,class(surfs{i}.vertices));

    vert_pos = 1;
    face_pos = 1;
    
    for i=1:n
        faces = size(surfs{i}.faces,1);
        verts = size(surfs{i}.vertices,1);
        
        surf.faces(face_pos:face_pos+faces-1,:) = surfs{i}.faces +...
            vert_pos - 1;
        surf.vertices(vert_pos:vert_pos+verts-1,:) = surfs{i}.vertices;
        
        face_pos = face_pos + faces;
        vert_pos = vert_pos + verts;
    end
end

