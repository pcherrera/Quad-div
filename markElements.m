function varargout = markElements(theta,varargin)

% markElements uses the (generalized) Doerfler marking criterion to mark
% elements for refinement.
%
% Usage: [MARK1 [,MARK2,...]] = markElements(THETA,INDICATOR1 [,INDICATOR2,...])
%    or  [MARK1 [,MARK2,...]] = markElements(THETA,RHO,INDICATOR1 [,INDICATOR2,...])
%
% Description: This file is part of the HILBERT program package for the
%              numerical solution of the Laplace equation with mixed
%              boundary conditions by use of BEM in 2D.
%               
%              Suppose that INDICATORX denotes a vector of refinement
%              indicators. Let INDICATORS denote the union of all refinement
%              indicators. The marking strategy introduced by Doerfler 1996 aims
%              to find a subset MARKED_INDICATORS of the entire set INDICATORS
%              such that 
%
%                THETA*sum(INDICATORS) <= sum(MARKED_INDICATORS)
%
%              The marking parameter THETA is strictly between 0 and 1.
%              INDICATOR1, INDICATOR2,... are vectors of refinement indicators. 
%              MARKED1, MARKED2,... are column vectors of indices j such that, 
%              e.g., INDICATOR1(j) formally belongs to the set
%              MARKED_INDICATORS. The optional parameter RHO between 0 and 1 
%              gives the percent of INDICATORS which is at least marked. For
%              RHO = 0 or RHO not given, the constructed MARKED_INDICATORS
%              has minimal cardinality.
%
%              Example: MARK = markElements(THETA,INDICATORS) provides the
%              minimal set of indices such that the Doerfler criterion is
%              satisfied.
%
%              Example: MARK = markElements(THETA,0.25,INDICATORS) provides a
%              set MARK that satisfies the Doerfler criterion and that contains
%              at least the 25% largest entries of INDICATORS.
%
% (C) 2009-2012 HILBERT-Team '09, '10, '12
% license + details:     http://www.asc.tuwien.ac.at/abem/hilbert
% support + bug report:  hilbert@asc.tuwien.ac.at
%
% Version: 3.0

%*** check whether optional parameter rho is given or not
if nargin == nargout +1
    rho = 0;
else
    rho = varargin{1};
    varargin = varargin(2:end);
end

%*** enforce input parameters to be column vectors and count their length
nE = zeros(1,nargout+1);
for j = 1:nargout
    nE(j+1) = length(varargin{j});
    varargin{j} = reshape(varargin{j},nE(j+1),1);
end

%*** generate set of all indicators
indicators = cat(1,varargin{:});
nE = cumsum(nE);

%*** realization of Doerfler marking
[indicators,idx] = sort(indicators,'descend');
sum_indicators = cumsum(indicators);
ell = max(ceil(rho*nE(end)),find(sum_indicators>=sum_indicators(end)*theta,1));
marked = idx(1:ell);

%*** split subset marked into subsets with respect to input vectors
for j = 1:nargout
    varargout{j} = marked( marked>nE(j) & marked<=nE(j+1) ) - nE(j);
end

