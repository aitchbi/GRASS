function Ap = ml_prune_adjacency(A,mask,surf,varargin)
% ML_PRUNE_ADJACENCY Remove edges from adjacency matrix using surface.
%   Ap = ML_PRUNE_ADJACENCY(A,mask,surf)
%       Prunes edges in the adjacency matrix A that intersect the provided
%       mesh surf. mask is a 3D tensor whose non-zero elements, indexed
%       linearly, correspond to the nodes in A.
%   Ap = ML_PRUNE_ADJACENCY(A,mask,surf,...)
%       Additional arguments can be provided as key-value pairs.
%         - conn6 (true/false) specifies whether or not A only has
%           6-connectivity. If 6-connectivity is used this can be set to
%           true for increased speed. For higher connectivity this must be
%           set to false (default).
%         - parallelize (true/false) specifies whether or not to run
%           prunings using multiple surfaces in parallel. Default is false.
%         - progess (true/false) specifies whether or not to report
%           progress during pruning. Default is false. If parallelization
%           is used no progress is reported.
%
%   Author:
%       Martin Larsson
%       March 2017

    control_params = {
        'conn6',0,...
        'parallelize',0,...
        'progress',0};

    conn6 = 0;
    parallelize = 0;
    progress = 0;

    argselectCheck(control_params,varargin);
    argselectAssign(varargin);

    if ~parallelize
        Ap = ml_prune_adjacency_private(A,mask,surf,conn6,progress);
    else
        % Divide faces evenly between workers.
        p = gcp;
        surfs = cell(p.NumWorkers,1);
        bounds = floor((0:p.NumWorkers)*size(surf.faces,1)/p.NumWorkers)+1;
        bounds(end) = size(surf.faces,1);

        for i=1:p.NumWorkers
            surfs{i}.faces = surf.faces(bounds(i):bounds(i+1),:);
            surfs{i}.vertices = surf.vertices;
        end

        % Parallel execution.
        spmd
            Aps = ml_prune_adjacency_private(A,mask,surfs{labindex},...
                conn6,false);
        end

        % Combine results from workers.
        Ap = Aps{1};
        for i=2:p.NumWorkers
            Ap = Ap & Aps{i};
        end
        Ap = Ap.*A;
    end
end

function Ap = ml_prune_adjacency_private(A,mask,g,conn6,progress)
% ML_PRUNE_ADJACENCY Remove edges from adjacency matrix using surface.
%   Ap = ML_PRUNE_ADJACENCY(A,mask,g)
%       Prunes edges in the adjacency matrix A that intersect the surface
%       g. mask is a volume indicating where the nodes of the graph are in
%       3D space. The non-zero elements in mask, indexed linearly,
%       correspond to the row/column numbers of A. g is a triangulated mesh
%       struct with the fields:
%           - vertices (Nx3)
%           - faces (Mx3)
%   Ap = ML_PRUNE_ADJACENCY(A,mask,g,conn6)
%       conn6 specifies whether or not A only has 6-connectivity. If
%       6-connectivity is used this can be set to true for increased speed.
%       For higher connectivity this must be set to false (default).
%   Ap = ML_PRUNE_ADJACENCY(A,mask,g,conn6,progress)
%       progress specifies whether or not to report progress. Default value
%       is false.
%
%   Author:
%       Martin Larsson
%       March 2017

    if nargin < 4
        conn6 = false;
    end
    if nargin < 5
        progress = false;
    end

    Ap = A;
    dim = size(mask);
    face_count = size(g.faces,1);
    
    indices = find(mask);
    indices_inv = zeros(dim);
    indices_inv(indices) = 1:length(indices);
    [indx,indy,indz] = ind2sub(dim,indices);
    
    % The size of to_remove might be a bit excessively large and arbitrary.
    % Not preallocating does not seem to significantly affect performance,
    % however it is kept in case the behavior is different for other
    % versions of MATLAB.
    to_remove = zeros(nnz(Ap),1);
    to_remove_index = 1;
    
    % It is faster to precalculate the triangles here than to extract the
    % relevant vertices in the for loop.
    all_tri = zeros(3,3,face_count);
    all_tri(1,:,:) = g.vertices(g.faces(:,1),:)';
    all_tri(2,:,:) = g.vertices(g.faces(:,2),:)';
    all_tri(3,:,:) = g.vertices(g.faces(:,3),:)';
    
    tri.faces = [1 2 3];

    if progress
        fprintf('* Pruning adjacency matrix...\nProgress: ');
        fprintf('%6d/%-6d\n',1,face_count);
    end
    
    % Handle one triangle at a time.
    for i = 1:face_count
        if progress && (~rem(i,500) || i == face_count)
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b%6d/%-6d\n',i,...
                face_count);
        end

        tri.vertices = all_tri(:,:,i);
        
        % Construct bounding box for the triangle.
        minv = floor(min(tri.vertices));
        maxv = ceil(max(tri.vertices));
        
        % Expand to accommodate dilation below.
        if ~conn6
            minv = minv - 1;
            maxv = maxv + 1;
        end
        
        minv = max(minv, 1);
        maxv = min(maxv, dim);
        
        bbdim = maxv - minv + 1;
        
        % Translate vertices from the mask space to the bounding box.
        tri.vertices = tri.vertices - repmat(minv,3,1) + 1;
        
        % Voxelization of triangle.
        V = polygon2voxel(tri, bbdim, 'none', false);

        % Perform dilation with 6-connectivity. This is not pretty but it
        % is apparently faster than using the imdilate function. This is
        % not necessary if A only has 6-connectivity, but is required to
        % catch "diagonal" edges.
        if ~conn6
            V2 = V;
            V2(1:end-1,:,:) = V2(1:end-1,:,:) | V(2:end,:,:);
            V2(:,1:end-1,:) = V2(:,1:end-1,:) | V(:,2:end,:);
            V2(:,:,1:end-1) = V2(:,:,1:end-1) | V(:,:,2:end);
            V2(2:end,:,:) = V2(2:end,:,:) | V(1:end-1,:,:);
            V2(:,2:end,:) = V2(:,2:end,:) | V(:,1:end-1,:);
            V2(:,:,2:end) = V2(:,:,2:end) | V(:,:,1:end-1);
            V = V2;
        end
        
        % Remove voxels not in the mask.
        V = V & mask(minv(1):maxv(1),minv(2):maxv(2),minv(3):maxv(3));
        
        % Translate voxels back to the mask space.
        [vx,vy,vz] = ind2sub(bbdim,find(V));
        if size(vx,1) == 0
            continue;
        end
        voxels = [vx vy vz];
        voxels = voxels + repmat(minv,size(voxels,1),1) - 1;
        
        % Convert mask space (x,y,z) to index in A.
        Aind = indices_inv(sub2ind(dim,voxels(:,1),voxels(:,2),...
            voxels(:,3)));
        
        % Define start and end point for all edges from the voxels.
        edge_index = 0;
        for j = 1:length(Aind)
            neighbors = find(Ap(:,Aind(j)));
            nn = length(neighbors);
            startv = repmat(voxels(j,:),nn,1);
            endv = [indx(neighbors) indy(neighbors) indz(neighbors)];
            
            orig(edge_index+1:edge_index+nn,:) = startv;
            dir(edge_index+1:edge_index+nn,:) = endv - startv;
            
            edges(edge_index+1:edge_index+nn) = sub2ind(size(Ap),...
                neighbors,repmat(Aind(j),nn,1));
            edge_index = edge_index + nn;
        end
        
        % There were no edges to check intersection with.
        if edge_index == 0
            continue;
        end

        % Check for intersections.
        intersections = TriangleRayIntersection(orig(1:edge_index,:),...
            dir(1:edge_index,:), all_tri(1,:,i), all_tri(2,:,i),...
            all_tri(3,:,i), 'lineType', 'segment', 'border', 'inclusive');
        
        edges = edges(intersections);
        
        % Append edges to list to be removed in bulk later.
        to_remove(to_remove_index+(0:length(edges)-1)) = edges;
        to_remove_index = to_remove_index + length(edges);
    end
    
    % Remove edges and enforce symmetry.
    Ap(to_remove(1:to_remove_index-1)) = 0;
    Ap = (Ap & Ap').*Ap;
end

