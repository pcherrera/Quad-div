function [coordinates,elements,boundary,edges,element2edges,...
    boundary2edges] = prepareFBMesh(coordinates,elements,boundary)

%*** sort Mesh sucht that boundary nodes are last!
[coordinates,elements,boundary] = sortMeshFEMBEM(coordinates,elements,...
    boundary);

%*** provide geometric data, i.e. edges, element2edges etc.
[edges, element2edges, boundary2edges] = ...
    provideGeometricData(elements,boundary);

%*** sort edges sucht that boundary elements, are last in the edges vector
[edges,element2edges,boundary2edges] = sortMeshFEMBEM(edges,...
    element2edges,boundary2edges);

%*** re-sort edges such that boundary2edges = (nEd-nB+1:nEd)'
%*** where nEd is number of edges and nB is number of boundary elements
nEd = size(edges,1);
nB = size(boundary,1);
[~,permutation] = sort(boundary2edges);

boundary2edges = boundary2edges(permutation);

permutation = [(1:nEd-nB)'; nEd-nB+permutation];

element2edges = permutation(element2edges);
edges = edges(permutation,:);

%*** make sure that edges(boundary2edges,:) == boundary
%*** this makes sense, so that edges have the right orientation on the
%boundary!
edges(boundary2edges,:) = boundary;
