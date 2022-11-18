function [edge2nodes,element2edges,varargout] ...
             = provideGeometricData(elements,varargin)

%provideGeometricData: returns geometric data for finite element mesh
%
%Usage:
%
%[EDGE2NODES,ELEMENT2EDGES,DIRICHLET2EDGES,NEUMANN2EDGES] ...
%    = PROVIDEGEOMETRICDATA(ELEMENTS,DIRICHLET,NEUMANN)
%
%Comments:
%
%    provideGeometricData expects as input a finite element mesh described
%    by the fields COORDINATES, ELEMENTS, DIRICHLET, and NEUMANN. The
%    function chooses a numbering of the edges and then return this numbering
%    related to nodes, edges, and the boundary conditions.
%
%    EDGE2NODES(K) returns the indices of the two nodes of the k-th edge.
%    ELEMENT2EDGES(J,K) provides the edge number of the edge between the two 
%    nodes ELEMENTS(J,K) and ELEMENTS(J,K+1). DIRICHLET2EDGES(K) provides
%    the number of the k-th Dirichlet edge given by DIRICHLET(K,:). The same
%    applies for NEUMANN2EDGES
%
%Remark:
%
%    This program is a supplement to the paper 
%    >> Efficient Implementation of Adaptive P1-FEM in Matlab <<
%    by S. Funken, D. Praetorius, and P. Wissgott. The reader should 
%    consult that paper for more information.   
%
%Authors:
% 
%    S. Funken, D. Praetorius, P. Wissgott  10-07-08

nE = size(elements,1);
nB = nargin-1;
%*** Node vectors of all edges (interior edges appear twice) 
I = elements(:);
J = reshape(elements(:,[2,3,1]),3*nE,1);
%*** Symmetrize I and J (so far boundary edges appear only once)
pointer = [1,3*nE,zeros(1,nB)];
for j = 1:nB
    boundary = varargin{j};
    if ~isempty(boundary)
        I = [I;boundary(:,2)];
        J = [J;boundary(:,1)];
    end
    pointer(j+2) = pointer(j+1) + size(boundary,1);
end
%*** Create numbering of edges
idxIJ = find(I < J);
edgeNumber = zeros(length(I),1);
edgeNumber(idxIJ) = 1:length(idxIJ);
idxJI = find(I > J);
number2edges = sparse(I(idxIJ),J(idxIJ),1:length(idxIJ));
[foo{1:2},numberingIJ] = find( number2edges );
[foo{1:2},idxJI2IJ] = find( sparse(J(idxJI),I(idxJI),idxJI) );
edgeNumber(idxJI2IJ) = numberingIJ;
%*** Provide element2edges and edge2nodes
element2edges = reshape(edgeNumber(1:3*nE),nE,3);
edge2nodes = [I(idxIJ),J(idxIJ)];
%*** Provide boundary2edges
for j = 1:nB
    varargout{j} = edgeNumber(pointer(j+1)+1:pointer(j+2));
end
