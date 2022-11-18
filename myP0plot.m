function myP0plot(coordinates,elements,x,varargin)

nE = size(elements,1);

if(nargin==4)
    p = varargin{1};
else
    p = 0;
end

if(p==0)
    % positions of all "Corners"
    C = coordinates(elements(:),:);
    % Scalar field defined at corners
    SC = repmat(x,size(elements,2),1);
    % new disjoint triangle indices into C for each original triangle
    FC = bsxfun(@plus,size(elements,1)*(0:size(elements,2)-1),(1:size(elements,1))');
    % plot each triangle as piecewise linear function (each of which happens to be constant)
    trisurf(FC,C(:,1),C(:,2),SC);
elseif(p==1)
    pP = (p+1)*(p+2)/2;
    % positions of all "Corners"
    C = coordinates(elements(:),:);
    % Scalar field defined at corners
    %SC = repmat(x,size(elements,2),1);
    SC = [x(1:pP:nE*pP); ...
          x(1:pP:nE*pP)+x(2:pP:nE*pP); ...
          x(1:pP:nE*pP)+x(3:pP:nE*pP) ];
    % new disjoint triangle indices into C for each original triangle
    FC = bsxfun(@plus,size(elements,1)*(0:size(elements,2)-1),(1:size(elements,1))');
    % plot each triangle as piecewise linear function (each of which happens to be constant)
    trisurf(FC,C(:,1),C(:,2),SC);
else
    error('pol. degree not supported');
end
    