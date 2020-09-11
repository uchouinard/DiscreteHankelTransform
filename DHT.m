function [F, Fgrid] = DHT(n, data, q, N)
% Calculate the Discrete Hankel Transform 
%
% This code is based on: Baddour, N. and Chouinard, U., 2017. Journal
% of Open Research Software, 5(1), p.4. doi.org/10.5334/jors.82
% with a few improvments:
% 1. The Y matrix code is now vectorized making it ~ x20 faster. 
% 2. Optional zero padding input similar to Matlab's fft functionality.
% 3. Supports also arrays similar to Matlab's fft functionality, meaning
%    that the Bessel zeros and Y matrix is only calculated once per run.
% 4. Updated bessel zeros calculation code (by Jason Nicholson).
%  
%
% Inputs:
%  n      - the order of the Hankel Transform
%  data   - the input data, a vector of size (q x 1), or a 2D array (q x n)
%  q      - the grid points of vec
%  N      - If the length of data is less than N, then data is padded with
%            trailing zeros to length N.
%
% Outputs:
%
%  F      - The transfomred vector
%  Fgrid  - The transformed grid
%
%   Ver 1 (2020-09-07)
%   Adi Natan (natan@stanford.edu)

%% example: generate array of signals in q range, pad zeros to N
%  N    = 2^8;               % zero padding
%  q    = eps:0.075:10;      % data grid
%  a    = linspace(1,20);     
%  data = sin(q'*a)./(q'*a); % set of data vectors
%  filt = sin(pi*[1:numel(q)]'./(numel(q))).^2;
%  data = bsxfun(@times,data, filt); % filter edges with sin^2 window
%  [F, Fgrid] = DHT(0.5, data, q, N);
% 
%  subplot(2,1,1)
%  imagesc(a,q,data);
%  xlim([a(1) a(end)]);
%  ylim([q(1) q(end)]);
%  title('data')
% 
%  subplot(2,1,2)
%  imagesc(a,Fgrid,F);
%  xlim([a(1) a(end)]);
%  ylim([a(1) a(end)]);
%  title('DHT(data)')


%% copyright of the original code in doi.org/10.5334/jors.82 that this code is based on:   
%{
_______________________________________________________________________
Copyright (c) 2015, Ugo Chouinard
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
________________________________________________________________________
%}

%% 
% Defaults
% if there is no padding argument
if (nargin < 4);  N=numel(q); end

%define grid resolution
dq = q(2) - q(1);

%create a vector of bessel zeroes:
BZ = bz(n,N);

%create the transformation matrix (yMatrix) using BZ
JN1 = besselj(n+1, BZ(1:N-1));    
M = 2*besselj(n,BZ(1:N-1)'*BZ(1:N-1)/BZ(N));
yMatrix = bsxfun(@rdivide, M, BZ(N)*JN1.^2');

% data grid given padding (if no padding qp=q).
qp = q(1):dq:(dq*(N-1));

% resample data on the grid defined by bessel zeros using resample_grid
resample_grid = @(R,Z) (Z(1:end-1)*R)/Z(end);

% resampled grid:
qs = resample_grid(qp(end), BZ);

% resample data according to qs:
data_s = zeros(numel(qs),size(data,2));
qrange = qs>=q(1) & qs<=q(end);
data_s(qrange,:) = interp1(q,data,qs(qrange),'spline');

%there is also the need to apply a scaling factor defined as:
scalingFactor = qp(end)^2/BZ(N); % q limited function scaling factor
 
% transform data (each data vector needs to be a column vector):
F = zeros(size(data_s));
for ii = 1:size(data_s,2)
    F(:,ii) = yMatrix*data_s(:,ii)*scalingFactor;
end

% the corresponding frequency grid is resampled with BZ:
resample_freq_grid = @(R,Z) Z(1:end-1)/R;

% resampled freq grid:
Fgrid = resample_freq_grid(qp(end),BZ);

end


function x = bz(n, k, kind)
%  calculates the zeros of Bessel function of the first and second kind
%
%   x = bz(n)
%   x = bz(n, k)
%   x = bz(n, k, kind)
%
%% Inputs
% * *n* - The order of the bessel function. n can be a scalar, vector, or
%       matrix.  n can be positive, negative, fractional, or any
%       combinaiton. abs(n) must be less than or equal to
%       146222.16674537213 or 370030.762407380 for first and second kind
%       respectively. Above these values, this algorithm will not find the
%       correct zeros because of the starting values therefore an error is
%       thrown instead.
% * k - The number of postive zeros to calculate.  When k is not supplied,
%       k = 5 is the default. k must be a scalar.
% * kind - kind is either 1 or 2. When kind is not supplied, default is
%          kind = 1.
%
%% Outputs
% * x - The calculated zeros.  size(x) = [size(n) k].
%
%% Description
% besselzero calculates the first k positive zeros of nth order bessel
% function of the first or second kind.  Note, that zero is not included as
% the first zero.
%
%% Algorithm
% the first three roots of any order bessel can be approximated by a simple
% equations.  These equations were generated using a least squares fit of
% the roots from orders of n=0:10000. The approximation is used to start
% the iteration of Halley's method.  The 4th and higher roots can be
% approximated by understanding the roots are regularly spaced for a given
% order.  Once the 2nd and 3rd roots are found, the spacing can be
% approximated by the distance between the 2nd and 3rd root.  Then again
% Halley's method can be applied to precisely locate the root.
%%
% Because the algorithm depends on good guesses of the first three zeros,
% if the guess is to far away then Halley's method will converge to the
% wrong zero which will subsequently cause any other zero to be incorrectly
% located. Therefore, a limit is put on abs(n) of 146222.16674537213 and
% 370030.762407380 for first and second kind respectively.  If n is
% specified above these limits, then an error is thrown.
%
%% Example
%   n = (1:2)';
%   k = 10;
%   kind = 1;
%   z = bz(n, k, kind);
%   x = linspace(0, z(end), 1000);
%   y = nan(2, length(x));
%   y(1,:) = besselj(n(1), x);
%   y(2,:) = besselj(n(2), x);
%   nz = nan(size(z));
%   nz(1,:) = besselj(n(1), z(1,:));
%   nz(2,:) = besselj(n(2), z(2,:));
%   plot(x, y, z, nz,'kx')
%
% Originally written by
% Written by: Greg von Winckel - 01/25/05
% Contact: gregvw(at)chtm(dot)unm(dot)edu
%
% Modified, Improved, and Documented by
% Jason Nicholson 2014-Nov-06
% Contact: jashale@yahoo.com
%% Change Log
% * Original release. 2005-Jan-25, Greg von Winckel.
% * Updated Documentation and commented algorithm. Fixed bug in finding the
%   the first zero of the bessel function of the second kind. Improved speed by
%   factor of 20. 2014-Nov-06, Jason Nicholson.
%
%% Input checking
assert(nargin >=1 | nargin <=3,'Wrong number of input arguments.');
% Take care of default cases of k and kind
if nargin < 2
    k = 5;
end
if nargin < 3
    kind = 1;
end
assert(isscalar(kind) & any(kind == [1 2]), '''kind''must be a scalar with value 1 or 2 only.');
assert(isscalar(k) & fix(k)==k & k>0, 'k must a positive scalar integer.');
assert(all(isreal(n(:))), 'n must be a real number.');
% negative orders have the same roots as the positive orders
n = abs(n);
% Check for that n is less than the ORDER_MAX
if kind==1
    ORDER_MAX = 146222.16674537213;
    assert(all(n <= ORDER_MAX), 'all n values must be less than or equal %6.10f for kind=1.', ORDER_MAX);
elseif kind==2
    ORDER_MAX = 370030.762407380;
    assert(all(n(:) <= ORDER_MAX), 'all n values must be less than or equal %6.10f for kind=2.', ORDER_MAX);
end
%% Setup Arrays
% output size
nSize = size(n);
if nSize(end) ==1
    outputSize = [nSize(1:end-1) k];
else
    outputSize = [nSize k];
end
% number of orders for each kth root
nOrdersPerRoot = prod(outputSize(1:end-1));
x = nan(outputSize);
%% Solve for Roots
switch kind
    case 1
        % coefficients and exponent are from least squares fitting the k=1,
        % n=0:10000.
        coefficients1j = [0.411557013144507;0.999986723293410;0.698028985524484;1.06977507291468];
        exponent1j = [0.335300369843979,0.339671493811664];
        % guess for k = 1
        x((1:nOrdersPerRoot)') = coefficients1j(1) + coefficients1j(2)*n(:) + coefficients1j(3)*(n(:)+1).^(exponent1j(1)) + coefficients1j(4)*(n(:)+1).^(exponent1j(2));
        % find first zero
        x((1:nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, 1, x0, kind), n(:), x((1:nOrdersPerRoot)'));
        
        if k >= 2
            % coefficients and exponent are from least squares fitting the k=2,
            % n=0:10000.
            coefficients2j = [1.93395115137444;1.00007656297072;-0.805720018377132;3.38764629174694];
            exponent2j = [0.456215294517928,0.388380341189200];
            % guess for k = 2
            x((nOrdersPerRoot+1:2*nOrdersPerRoot)') = coefficients2j(1) + coefficients2j(2)*n(:) + coefficients2j(3)*(n(:)+1).^(exponent2j(1)) + coefficients2j(4)*(n(:)+1).^(exponent2j(2));
            % find second zero
            x((nOrdersPerRoot+1:2*nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, 2, x0, kind), n(:), x((nOrdersPerRoot+1:2*nOrdersPerRoot)'));
        end
        
        if k >= 3
            % coefficients and exponent are from least squares fitting the k=3,
            % n=0:10000.
            coefficients3j = [5.40770803992613;1.00093850589418;2.66926179799040;-0.174925559314932];
            exponent3j = [0.429702214054531,0.633480051735955];
            % guess for k = 3
            x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)') = coefficients3j(1) + coefficients3j(2)*n(:) + coefficients3j(3)*(n(:)+1).^(exponent3j(1)) + coefficients3j(4)*(n(:)+1).^(exponent3j(2));
            % find second zero
            x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, 3, x0, kind), n(:), x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)'));
        end
    case 2
        % coefficients and exponent are from least squares fitting the k=1,
        % n=0:10000.
        coefficients1y = [0.0795046982450635;0.999998378297752;0.890380645613825;0.0270604048106402];
        exponent1y = [0.335377217953294,0.308720059086699];
        % guess for k = 1
        x((1:nOrdersPerRoot)') = coefficients1y(1) + coefficients1y(2)*n(:) + coefficients1y(3)*(n(:)+1).^(exponent1y(1)) + coefficients1y(4)*(n(:)+1).^(exponent1y(2));
        % find first zero
        x((1:nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, 1, x0, kind), n(:), x((1:nOrdersPerRoot)'));
        
        if k >= 2
            % coefficients and exponent are from least squares fitting the k=2,
            % n=0:10000.
            coefficients2y = [1.04502538172394;1.00002054874161;-0.437921325402985;2.70113114990400];
            exponent2y = [0.434823025111322,0.366245194174671];
            % guess for k = 2
            x((nOrdersPerRoot+1:2*nOrdersPerRoot)') = coefficients2y(1) + coefficients2y(2)*n(:) + coefficients2y(3)*(n(:)+1).^(exponent2y(1)) + coefficients2y(4)*(n(:)+1).^(exponent2y(2));
            % find second zero
            x((nOrdersPerRoot+1:2*nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, 2, x0, kind), n(:), x((nOrdersPerRoot+1:2*nOrdersPerRoot)'));
        end
        
        if k >= 3
            % coefficients and exponent are from least squares fitting the k=3,
            % n=0:10000.
            coefficients3y = [3.72777931751914;1.00035294977757;2.68566718444899;-0.112980454967090];
            exponent3y = [0.398247585896959,0.604770035236606];
            % guess for k = 3
            x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)') = coefficients3y(1) + coefficients3y(2)*n(:) + coefficients3y(3)*(n(:)+1).^(exponent3y(1)) + coefficients3y(4)*(n(:)+1).^(exponent3y(2));
            % find second zero
            x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, 3, x0, kind), n(:), x((2*nOrdersPerRoot+1:3*nOrdersPerRoot)'));
        end
    otherwise
        error('Code should never get here.');
end
if k >= 4
    for iRoot = 4:k
        % guesses for remaining roots x[k] = rootSpacing + x[k-1]
        x(((iRoot-1)*nOrdersPerRoot+1:iRoot*nOrdersPerRoot)') = 2*x(((iRoot-2)*nOrdersPerRoot+1:(iRoot-1)*nOrdersPerRoot)')- x(((iRoot-3)*nOrdersPerRoot+1:(iRoot-2)*nOrdersPerRoot)');
        % find the remaining zeros
        x(((iRoot-1)*nOrdersPerRoot+1:iRoot*nOrdersPerRoot)') = arrayfun(@(n, x0) findzero(n, k, x0, kind), n, x(((iRoot-1)*nOrdersPerRoot+1:iRoot*nOrdersPerRoot)'));
    end
end
end
function x=findzero(n,k,x0,kind)
% Uses Halley's method to find a zero given the starting point x0
% http://en.wikipedia.org/wiki/Halley's_method
ITERATIONS_MAX = 100;       % Maximum number of iteration
TOLERANCE_RELATIVE = 1e4;   % 16-4 = 12 significant digits
% Setup loop
error = 1;
loopCount = 0;
x = 1; % Initialization value only.  It is is not used.
% Begin loop
while abs(error)>eps(x)*TOLERANCE_RELATIVE && loopCount<ITERATIONS_MAX
    
    switch kind
        case 1
            a = besselj(n,x0);
            b = besselj((n+1),x0);
        case 2
            a = bessely(n,x0);
            b = bessely((n+1),x0);
    end
    
    xSquared = x0*x0;
    
    error = 2*a*x0*(n*a-b*x0)/(2*b*b*xSquared-a*b*x0*(4*n+1)+(n*(n+1)+xSquared)*a*a);
    
    % Prepare for next loop
    x=x0-error;
    x0=x;
    loopCount=loopCount+1;
    
end
% Handle maximum iterations
if loopCount>ITERATIONS_MAX-1
    warning('Failed to converge to within relative tolerance of %e for n=%f and k=%d in %d iterations', eps(x)*TOLERANCE_RELATIVE, n, k, ITERATIONS_MAX);
end
end
