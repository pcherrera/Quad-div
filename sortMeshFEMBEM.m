function [coordinates,elements,boundary,varargout] = ... 
    sortMeshFEMBEM(coordinates,elements,boundary,varargin)
%SORTMESHFEMBEM    sorts Finite Element Mesh such that:
%
%   -) boundary nodes are the last nodes in coordinates
% 
%   
% TF: 23/04/2013

nE = size(elements,1);

%*** determine nodes on boundary (Dirichlet) and interior nodes (freenodes)
nC = size(coordinates,1);
nodes_boundary = unique(boundary);
freenodes = setdiff((1:nC)',nodes_boundary);

%*** build permutation such that Dirichlet nodes are last
nodes = [freenodes;nodes_boundary];
[foo,permutation] = sort(nodes);

%*** permute indices of nodes
coordinates(permutation,:) = coordinates;
elements = permutation(elements);
boundary = permutation(boundary);

if(length(varargin)==1)
    Tmat = varargin{1};
    Tmat(:,permutation) = Tmat;
    varargout{1} = Tmat;
end