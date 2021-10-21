function Vt = ml_transform_vertices(A,V)
% ML_TRANFORM_VERTICES Transform 3D vertices.
%   Vt = ML_TRANFORM_VERTICES(A,V)
%       Transforms the 3D vertices in V using the 4x4 matrix A. V can be
%       either 3xN or Nx3. Vt is the transformed vertices and has the same
%       size as V.
%
%   Author:
%       Martin Larsson
%       March 2017

    if size(V,1) == 3
        Vt = A * [V; ones(1,size(V,2))];
        Vt = Vt(1:3,:);
    elseif size(V,2) == 3
        Vt = (A * [V'; ones(1,size(V,1))])';
        Vt = Vt(:,1:3);
    else
        error('V must have either 3 rows or 3 columns.');
    end
end

