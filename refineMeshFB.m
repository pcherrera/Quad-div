function [vertices,new_volumes,varargout] = ...
             refineMesh(vertices,volumes,varargin)
% TF: ADAPTION FOR MULTILEVEL ALGORITHMS
% PROVIDES Transfermatrix for old2newcoordinates
% 10/04/2013

% refineMesh refines a given volume mesh together with its restriction
% to the boundary via Newest Vertex Bisection (NVB)
% of certain marked volumes and boundary edges. Furthermore it
% guarantees that the resulting mesh is regular. Optional refineMesh
% transferMatates the coarse mesh solution xold on the fine mesh (nodal
% interpolation on the volumes and constant continuation on the boundary).
% 
% Usage: [VERTICES,VOLUMES,BOUNDARY [[,FATHER2VOLUMES,FATHER2BOUNDARY][,xnew]] ...
%            = refineMesh(VERTICES,VOLUMES,BOUNDARY ...
%                          [[,MARKED_VOLUMES,MARKED_BOUNDARY][,xold]])
%
%    or  [VERTICES,VOLUMES,DIRICHLET,NEUMANN ...
%               [[,FATHER2VOLUMES,FATHER2DIRICHLET,FATHER2NEUMANN][,xnew]]] ...
%            = refineMesh(VERTICES,VOLUMES,DIRICHLET,NEUMANN ...
%                          [,MARKED_VOLUMES,MARKED_DIRICHLET,MARKED_NEUMANN][,xold]])
%
%    or  [VERTICES,VOLUMES,<BOUNDARIES> ...
%               [[,FATHER2VOLUMES,<FATHER2BOUNDARIES>][,xnew]]] ...
%            = refineMesh(VERTICES,VOLUMES,<BOUNDARIES> ...
%                          [,MARKED_VOLUMES,<MARKED_BOUNDARIES>][,xold]])
%
%    where BOUNDARIES is a synonym for DIRICHLET, NEUMANN, ROBIN, ...
%    so that one can manage any finite number of boundary parts
%
% Description: This file is part of the HILBERT program package for the
%              numerical solution of the Laplace equation with mixed 
%              boundary conditions by use of BEM or FEM-BEM in 2D.
%
%              Let {T1,...,TN} be a mesh described by matrices VERTICES and
%              VOLUMES. The restriction to the boundary {E1,...,Ek} is 
%              described in terms of matrices BOUNDARIES which represent 
%              different boundary parts.
%              The optional (m x 1)-vector MARKED_VOLUMES 
%              contains the indices of all volumes Tj which will be refined by
%              newest vertex bisection NVB whereas the optional (ni x 1)-vectors
%              MARKED_BOUNDARIES hold the indices of all boundary elements Ek 
%              which will be refined by bisection. The refined triangulation is
%              returned in terms of extended matrices VERTICES, VOLUMES and
%              BOUNDARIES, where the latter ones hold again the restriction
%              to the boundary, i.e. the indices of vertices which represent 
%              the boundary elements. Moreover, the optional (N x 4)-matrix
%              FATHER2VOLUMES links the initial mesh with the refined mesh
%              in the sense that FATHER2VOLUMES(j,:) contains the indices of
%              the fine volumes which are the sons of the coarse-mesh volume Tj.
%              We distinguish different cases:
%              1. If Tj has not been refined FATHER2VOLUMES(j,:) contains
%              four times the index of Tj with respect to the refined mesh.
%              2. If Tj is refined by NVB of the first edge 
%              FATHER2VOLUME(j,:) contains the indices of its two sons, each 
%              of it appears twice, i.e. 
%              FATHER2VOLUMES(j,1)=FATHER2VOLUMES(j,2) and
%              FATHER2VOLUMES(j,3)=FATHER2VOLUMES(j,4).
%              3. If Tj is refined by NVB of the first and second edge
%              FATHER2VOLUME(j,:) contains the indices of its three sons, and 
%              there holds FATHER2VOLUMES(j,1)=FATHER2VOLUMES(j,2). 
%              4. If Tj is refined by NVB of the first and third edge
%              FATHER2VOLUME(j,:) contains the indices of its three sons, and 
%              there holds FATHER2VOLUMES(j,3)=FATHER2VOLUMES(j,4). 
%              5. If Tj is refined by NVB of all edges FATHER2VOLUME(j,:)
%              holds the indices of the four sons.
%              The optional (N x 2)-matrices FATHER2BOUNDARIES link the 
%              initial mesh elements of each boundary part to the refined mesh
%              elements in the sense that FATHER2BOUNDARIES(k,:)
%              contains the indices of the fine elements which are the sons of
%              the coarse-mesh element Ek. If Ek has not been refined,
%              FATHER2BOUNDARIES(j,1) = FATHER2BOUNDARIES(j,2) contains the index of 
%              Tj with respect to the refined mesh.
%              
%              If the optional vectors MARKED_VOLUMES and MARKED_BOUNDARIES
%              are not given, the function performs a uniform refinement,
%              i.e., all volumes Tj and all boundary elements Ek are refined.
%
%              xold is the solution on the coarse mesh, represented by an
%              array of the form [nC x 1].
%                                [nB x 1]
%              Where nC is the number of vertices and nB is the number
%              of boundary elements.                                            
%
% Author: Markus Aurada
%
% (C) M. Aurada, M. Ebner, S. Ferraz-Leite,
%     P. Goldenits, M. Karkulik, M. Mayr,
%     D. Praetorius
%
% Version: 0.5

%*** count number of boundary parts from input 
%*** nB will hold this number 
%*** nB_elts will hold the number of elements of each boundary part
nB = 0;

for iter = 1 : (nargin - 2)
    if size(varargin{iter},2) == 2
        nB = nB + 1;
        nB_elts(iter) = size(varargin{iter},1);
    else
        break;
    end

end
%*** check if transferMatation of the coarse solution is wanted
if((nB + 3==nargin) || (2*nB+4==nargin))
    transferMat=true;
else
    transferMat=false;
end
%*** check if there is at least one boundary part
if nB == 0
   error('refineMesh: There has to be at least one boundary part!');
end

%*** check the correct number of input parameters
if ~( (nargin == (nB+2)) || (nargin == (2*nB+3)) || (nargin == nB+3) || (nargin== 2*nB+4))
    error('refineMesh: Wrong number of input arguments!');
end

%*** check the correct number of output parameters
if ~( (nargout == (nB+2)) || (nargout == (2*nB+3)) || (nargout == nB+3) || (nargout== 2*nB+4))
    error('refineMesh: Wrong number of output arguments!');
end

%*** check, if user asks for father2son fields in output
if nargout == (2*nB+3+transferMat)
    output_father2son = true;
else
    output_father2son = false;
end

%*** check, if user asks for adaptive refinement
if nargin == (2*nB+3+transferMat) %***Achtung Änderung: vorher, if nargout == (2*nB+3)   
   adaptive = true;
else
   adaptive = false;
end

%*** number of volumes
nVols = size(volumes,1);

%*** Obtain geometric information on edges
%*** Node vectors of all edges (interior edges appear twice) 
I = volumes(:);
J = reshape(volumes(:,[2,3,1]),3*nVols,1);

%*** obtain set of all boundary elements of the boundary partition
bdry_elts = cat(1,varargin{1 : nB});

%*** indices of boundary parts w.r.t. entire boundary elements
%*** and w.r.t. entire edges
ptr_bdry_elts = cumsum([0,nB_elts]);
ptr_bdry_edges = ptr_bdry_elts + 3*nVols;

%*** Node vectors of all edges (all edges appear twice) 
I = [I;bdry_elts(:,2)];
J = [J;bdry_elts(:,1)];

%*** Create numbering of edges
idx_IJ = find(I < J);
edge_number = zeros(length(I),1);
edge_number(idx_IJ) = 1:length(idx_IJ);
idx_JI = find(I > J);
nodes2edges = sparse(I(idx_IJ),J(idx_IJ),1:length(idx_IJ));
[foo{1:2},numbering_IJ] = find(nodes2edges);
[foo{1:2},idx_JI2IJ] = find(sparse(J(idx_JI),I(idx_JI),idx_JI));
edge_number(idx_JI2IJ) = numbering_IJ;

%*** Provide bdry_edges
for j = 1:nB
    bdry_edges{j} = edge_number(ptr_bdry_edges(j)+1:ptr_bdry_edges(j+1));
end

%*** Provide volumes2edges and edge2nodes
volumes2edges = reshape(edge_number(1:3*nVols),nVols,3);
edge2nodes = [I(idx_IJ),J(idx_IJ)];

%*** 1. determine whether uniform or adaptive mesh-refinement
%*** 2. in case of adaptive mesh-refinement compute vector 
%***    of marked volumes and vector of marked boundary edges.
%***    The latter will be computed w.r.t. entire boundary edges
if adaptive
   marked_bdry_elts = zeros(0,1); 
   for iter = 1 : nB
       marked_bdry_elts = [marked_bdry_elts; varargin{iter + nB + 1} + ...
                           ptr_bdry_elts(iter)];
   end

   marked_vols = varargin{nB + 1};
else
   marked_bdry_elts = [1:ptr_bdry_elts(end)]';
   marked_vols = [1:nVols]';
end

%*** compute index vectors for marking of boundary edges
bdry_nodes_min = min(bdry_elts,[],2);
bdry_nodes_max = max(bdry_elts,[],2);

%*** Mark edges for refinement
edge2new_node = zeros(length(edge2nodes),1);
edge2new_node(nodes2edges(sub2ind(size(nodes2edges),...
   bdry_nodes_min(marked_bdry_elts),bdry_nodes_max(marked_bdry_elts)))) = 1;
edge2new_node(volumes2edges(marked_vols,:)) = 1; %****Achtung Änderung: vorher edge2new_node(volumes2edges(marked_vols,1)) = 1; 
swap = 1;
while ~isempty(swap)
    marked_edges = edge2new_node(volumes2edges);
    swap = find( ~marked_edges(:,1) & (marked_edges(:,2) | marked_edges(:,3)) );
    edge2new_node(volumes2edges(swap,1)) = 1;
end
nCold=size(vertices,1);
%*** Generate new nodes
idx = find(edge2new_node);
edge2new_node(idx) = size(vertices,1) + (1:nnz(edge2new_node));
nV_old = size(vertices,1); % TF:
vertices(edge2new_node(idx),:) = ...
   0.5*(vertices(edge2nodes(idx,1),:)+vertices(edge2nodes(idx,2),:));

% %*** Generate values for xnew
% xnew=zeros(2,1);
% if(transferMat)
%    xnew(edge2new_node(idx))=0.5*(xold(edge2nodes(idx,1))+xold(edge2nodes(idx,2))); 
%    xnew(1:nCold)=xold(1:nCold);
% end

if(transferMat)
    %% Create Transfermatrix
    % Benuetze: Alte Knoten kommen genau am Anfang des Arrays der neuen Knoten
    % vor (mit gleichen Indizes)

    idx_e2nn = find(edge2new_node);
    nNn = size(idx_e2nn,1);

    % I_Tmat = (1:nV_old)';
    % J_Tmat = (1:nV_old)';
    % Tmat = ones(nV_old,1);
    % %*** Durchlaufe alle alten Knoten und suche Kanten zu jedem alten Knoten
    % %*** die zugehoerigen Kanten!
    % for j=1:nNn
    %     J_Tmat(nV_old+j) = edge2new_node(idx_e2nn(j));
    %     J_Tmat(nV_old+nNn+j) = edge2new_node(idx_e2nn(j));
    %     I_Tmat(nV_old+j) = edge2nodes(idx_e2nn(j),1);
    %     I_Tmat(nV_old+nNn+j) = edge2nodes(idx_e2nn(j),2);
    %     
    %     Tmat(nV_old+j) = 0.5;
    %     Tmat(nV_old+nNn+j) = 0.5;
    % end
    % Tmat = sparse(I_Tmat,J_Tmat,Tmat,nV_old,nV_old+nNn);


    % vectorized
    I_Tmat = [(1:nV_old)'; edge2nodes(idx_e2nn,1); edge2nodes(idx_e2nn,2)];
    J_Tmat = [(1:nV_old)'; edge2new_node(idx_e2nn); edge2new_node(idx_e2nn)];
    Tmat = sparse(I_Tmat,J_Tmat,[ones(nV_old,1); 0.5*ones(2*nNn,1)]);
    varargout{nargout-2} = Tmat;
end
%*** Refine boundary edges and build father2boundaries if asked for
for j = 1:nB
    bdry = bdry_elts(ptr_bdry_elts(j)+1:ptr_bdry_elts(j+1),:);
    new_nodes = edge2new_node(bdry_edges{j});
    marked_edges = find(new_nodes);
    nMkd_edgs = length(marked_edges);
    non_marked_edges = find(~new_nodes);
    nNon_mkd_edgs = length(non_marked_edges);
    if ~isempty(marked_edges)
        bdry = [bdry(non_marked_edges,:); ...
                bdry(marked_edges,1),new_nodes(marked_edges); ...
                new_nodes(marked_edges),bdry(marked_edges,2)];
%         if(transferMat)
%             xnew=[xnew;xold(nCold+non_marked_edges);xold(nCold+marked_edges);xold(nCold+marked_edges)];
%         end
    end
    varargout{j} = bdry;
    if output_father2son
       father2son = zeros(nB_elts(j),2);
       father2son([non_marked_edges;marked_edges],:) = ...
          [ repmat([1:nNon_mkd_edgs]',1,2);...
            [(nNon_mkd_edgs+1):(nNon_mkd_edgs+nMkd_edgs)]',...
            [(nNon_mkd_edgs+nMkd_edgs+1):(nNon_mkd_edgs+2*nMkd_edgs)]' ];
       varargout{nB+j+1} = father2son;
    end
end

%*** Provide new nodes for refinement of volumes
new_nodes = edge2new_node(volumes2edges);

%*** Determine type of refinement for each volume
marked_edges = (new_nodes~=0);
none = ~marked_edges(:,1);
bisec1   = ( marked_edges(:,1) & ~marked_edges(:,2) & ~marked_edges(:,3) );
bisec12  = ( marked_edges(:,1) &  marked_edges(:,2) & ~marked_edges(:,3) );
bisec13  = ( marked_edges(:,1) & ~marked_edges(:,2) &  marked_edges(:,3) );
bisec123 = ( marked_edges(:,1) &  marked_edges(:,2) &  marked_edges(:,3) );

%*** Generate volume numbering for refined mesh
idx = ones(nVols,1);
idx(bisec1)   = 2; %*** bisec(1): newest vertex bisection of 1st edge
idx(bisec12)  = 3; %*** bisec(2): newest vertex bisection of 1st and 2nd edge
idx(bisec13)  = 3; %*** bisec(2): newest vertex bisection of 1st and 3rd edge
idx(bisec123) = 4; %*** bisec(3): newest vertex bisection of all edges

idx = [1;1+cumsum(idx)];

%*** Generate new elements
new_volumes = zeros(idx(end)-1,3);
new_volumes(idx(none),:) = volumes(none,:);
new_volumes([idx(bisec1),1+idx(bisec1)],:) = ...
   [volumes(bisec1,3),volumes(bisec1,1),new_nodes(bisec1,1); ...
    volumes(bisec1,2),volumes(bisec1,3),new_nodes(bisec1,1)];
new_volumes([idx(bisec12),1+idx(bisec12),2+idx(bisec12)],:) = ...
   [volumes(bisec12,3),volumes(bisec12,1),new_nodes(bisec12,1); ...
    new_nodes(bisec12,1),volumes(bisec12,2),new_nodes(bisec12,2); ...
    volumes(bisec12,3),new_nodes(bisec12,1),new_nodes(bisec12,2)]; 
new_volumes([idx(bisec13),1+idx(bisec13),2+idx(bisec13)],:) = ...
   [new_nodes(bisec13,1),volumes(bisec13,3),new_nodes(bisec13,3); ...
    volumes(bisec13,1),new_nodes(bisec13,1),new_nodes(bisec13,3); ...
    volumes(bisec13,2),volumes(bisec13,3),new_nodes(bisec13,1)];
new_volumes([idx(bisec123),1+idx(bisec123),2+idx(bisec123),...
             3+idx(bisec123)],:) = ...
   [new_nodes(bisec123,1),volumes(bisec123,3),new_nodes(bisec123,3); ...
    volumes(bisec123,1),new_nodes(bisec123,1),new_nodes(bisec123,3); ...
    new_nodes(bisec123,1),volumes(bisec123,2),new_nodes(bisec123,2); ...
    volumes(bisec123,3),new_nodes(bisec123,1),new_nodes(bisec123,2)];

%*** build father2volumes
if output_father2son
   father2son = zeros(nVols,4);
   father2son(none,:) = repmat(idx(none),1,4);
   father2son(bisec1,:) = [repmat(idx(bisec1),1,2),repmat(idx(bisec1)+1,1,2)];
   father2son(bisec12,:) = ...
      [repmat(idx(bisec12),1,2),idx(bisec12)+1,idx(bisec12)+2];
   father2son(bisec13,:) = ...
      [idx(bisec13),idx(bisec13)+1,repmat(idx(bisec13)+2,1,2)];
   father2son(bisec123,:) = ...
      [idx(bisec123),idx(bisec123)+1,idx(bisec123)+2,idx(bisec123)+3];
   varargout{nB+1} = father2son;
end

% %*** sorting vertices such that boundary nodes appear first
% bdry_nodes = zeros(0,1); 
% for j = 1:nB
%     bdry_nodes = [bdry_nodes;setdiff(unique(varargout{j}),bdry_nodes)];
% end
% 
% %*** number of boundary nodes and number of vertices
% nBdry_nodes = length(bdry_nodes);
% nVerts = size(vertices,1);
% 
% %*** renumeration of vertices 
% idx = [bdry_nodes;setdiff([1:nVerts]',bdry_nodes)]; % sort nodes first
% verts2new_verts(idx) = ...
%    [ [1:nBdry_nodes],[(nBdry_nodes+1):nVerts] ]; % sort first
% 
% %*** reorder vertices
% vertices = vertices(idx,:);
% if(transferMat)   
%     xnew(1:size(vertices,1))=xnew(idx);
%     varargout{nargout-2}=xnew; 
% end
% %%% assign new vertice indices to new_volumes and to all boundary parts
% new_volumes=verts2new_verts(new_volumes);
% for j = 1:nB
%     varargout{j} = verts2new_verts(varargout{j});
% end
