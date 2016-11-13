function [V,Vy,s] = pppiv(V,Vy,varargin)

%PPPIV Post-processing of 2-D PIV velocity field.
%   [Vx2,Vy2] = PPPIV(Vx1,Vy1) carries out robust post-processing of 2-D
%   PIV velocity data. Vx1 and Vy1 must be two matrices of same size that
%   contain the x- and y-components of the velocities at equally spaced
%   points in the Cartesian plane.
%
%   PPPIV uses a robust penalized least squares method (acronymed DCT-PLS)
%   that makes the smoothed output (Vx2,Vy2) no dependent upon the outlying
%   (spurious) vectors: outliers are replaced and velocity vectors are
%   smoothed using a single automated process.
%
%   MISSING DATA: Non finite (NaN or Inf) values in (Vx1,Vy1) are
%   considered as missing velocities. The algorithm replaces them
%   automatically. 
%
%   PPPIV uses SMOOTHN. See <a href="matlab:help smoothn">SMOOTHN</a> for more details.
%
%   [Vx2,Vy2] = PPPIV(Vx1,Vy1,ROI) post-processes (Vx1,Vy1) in the region
%   of interest defined by the binary matrix ROI: 1 => inside the region of
%   interest, 0 => outside (i.e. masked data). NaN values are assigned to
%   (Vx2,Vy2) velocities outside the region of interest.
%
%   By default, PPPIV selects the smoothing parameter automatically by
%   minimizing the GCV score (see reference #1 for details). Alternatively,
%   the amount of smoothing can be somewhat adjusted by adding one of the
%   three smoothing options:
%   [...] = PPPIV(Vx1,Vy1,OPTION) or [...] = PPPIV(Vx1,Vy1,ROI,OPTION)
%   The available options are:
%          '2x2'            - weak smoothing
%          '3x3'            - medium smoothing
%          'nosmoothing'    - extremely weak smoothing
%   The DCT-PLS behaves similarly to a 2x2 and 3x3 moving average (in terms
%   of cut-off frequency at -3 dB, see reference #2) with the '2x2' and
%   '3x3' options, respectively. If you prefer to process the outliers
%   without smoothing the data then choose 'nosmoothing'.
%
%   The smoothing parameter can also be tuned manually using the following
%   syntax:
%   [...] = PPPIV(Vx1,Vy1,S) or [...] = PPPIV(Vx1,Vy1,ROI,S)
%   where S must be a scalar >0.
%
%   [Vx2,Vy2,S] = PPPIV(...) also returns the smoothness parameter S.
%
%   PPPIV (no input/output argument) runs the following example.
%
%   Notes
%   -----
%   Several Matlab functions available <a
%   href="matlab:web('http://www.biomecardio.com/matlab')">here</a> or <a
%   href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/71589')">here</a> are required:
%   SMOOTHN, DCTN, IDCTN
%
%   References
%   ---------- 
%   1) Garcia D, Robust smoothing of gridded data in one and higher
%   dimensions with missing values. Computational Statistics & Data
%   Analysis, 2010;54:1167-1178. <a
%   href="matlab:web('http://www.biomecardio.com/pageshtm/publi/csda10.pdf')">Link to the paper</a>.
%   2) Garcia D, A fast all-in-one method for post-processing of PIV data.
%   Exp in Fluids, 2010; under review.
%
%   Example:
%   -------
%   % -- Cellular vortical flow --
%   [x,y] = meshgrid(linspace(0,1,64));
%   Vx = cos(2*pi*x+pi/2).*cos(2*pi*y);
%   Vy = sin(2*pi*x+pi/2).*sin(2*pi*y);
%   % -- Corrupt the original flow --
%   Vx = Vx + sqrt(0.1)*randn(size(Vx)); % adding Gaussian noise
%   Vy = Vy + sqrt(0.1)*randn(size(Vx)); % adding Gaussian noise
%   I = randperm(numel(Vx));
%   n = round(numel(Vx)/5);
%   Vx(I(1:n)) = (rand(n,1)-0.5)*4; % adding outliers
%   Vy(I(1:n)) = (rand(n,1)-0.5)*4; % adding outliers
%   % -- Smooth the corrupt velocity field using PPPIV --
%   [Vx2,Vy2] = pppiv(Vx,Vy);
%   % -- Display the results --
%   figure, quiver(x,y,Vx,Vy,2.5), axis square
%   title('Noisy velocity field')
%   figure, quiver(x,y,Vx2,Vy2,2.5), axis square
%   title('Post-processed velocity field')
%
%   See also SMOOTHN
%
%   -- Damien Garcia -- 2009/10, revised 2010/06
%   -- Website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>

% SMOOTHN is required: test if SMOOTHN exists
text = ['SMOOTHN is required. Download SMOOTHN <a href="matlab:web(''',...
    'http://www.biomecardio.com/matlab/smoothn.html'')">here</a>.'];
assert(logical(exist('smoothn','file')),text)

if nargin==0&&nargout==0, RunTheExample, return, end

% Check input arguments
narginchk(2,4);
sizV = size(V);
opt = '';
ROI = true(sizV);
s = [];
if nargin==3
    if isnumeric(varargin{1}) && isscalar(varargin{1}) % PPPIV(Vx,Vy,s)
        s = varargin{1};
        opt = 's';
    elseif isnumeric(varargin{1}) || islogical(varargin{1}) % PPPIV(Vx,Vy,ROI)
        ROI = logical(varargin{1});
    elseif ischar(varargin{1}) % PPPIV(Vx,Vy,option)
        opt = varargin{1};
    end
elseif nargin==4 % PPPIV(Vx,Vy,ROI,option)
    ROI = logical(varargin{1});
    opt = varargin{2};
end
assert(isequal(sizV,size(Vy),size(ROI)),...
    'Input arrays must have same size.')

% Smoothing parameter option
switch lower(opt)
    case ''
    case 's', s = abs(s)+eps;
    case {'2x2','22','2'}, s = 0.10;
    case {'3x3','33','3'}, s = 0.54;
    case {'nosmoothing','no','nosmooth','nosm'}, s = 0.001;
    otherwise
        warning('MATLAB:pppiv:BadOption',...
            'The smoothing option must be: ''2x2'', ''3x3'' or ''nosmoothing''.')
end

% Complex velocity field
V = complex(V,Vy);

% Missing & masked data
V(~isfinite(V)) = Inf;
V(~ROI) = NaN;

% Remove surrounding NaNs surrounding the ROI
[V,Inanc] = nancrop(V);

% Robust smoothing using the DCT-PLS
wstate1 = warning('query','MATLAB:smoothn:MaxIter');
wstate2 = warning('query','MATLAB:smoothn:SLowerBound');
warning('off','MATLAB:smoothn:MaxIter')
warning('off','MATLAB:smoothn:SLowerBound')
% ----
if s<0.1 % "no smoothing" option
    V = smoothn(V,0.1,'robust'); % Coarse smoothing for initial conditions
    V = smoothn(V,s,'robust','TolZ',1e-4,'MaxIter',100,'Initial',V); 
else
    [V,s] = smoothn(V,s,'robust');
end
% ----
warning(wstate1)
warning(wstate2)

% Reconstruct the output velocities
Vy = (NaN+1i*NaN)*Inanc;
Vy(Inanc) = V;
Vy(~ROI) = NaN+1i*NaN;
V = real(Vy);
Vy = imag(Vy);

% Warning message in case of potential oversmoothing
if s>0.55
    warning('MATLAB:pppiv:LowSNR',...
        ['The optimal smoothness parameter is higher than 0.55 (s = ',...
        num2str(s,'%.2e') '). The PIV data have very likely a low SNR.'])
end

end

%% NANCROP
function [A,I] = nancrop(A)

%NANCROP Discard NaN-space around arrays.
%   B = NANCROP(A) eliminates NaN-space surrounding the array A. A can be
%   a vector, a matrix or a N-dimensional array.
%
%   [B,I] = NANCROP(A) also returns a boolean array I of same size as A
%   which represents the area taken up by B and so that B(:) = A(I)
%
%   -- Damien Garcia -- 2008/07, revised 2009/11
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>

I = true(size(A));
J = isnan(A);

for i = 1:ndims(A)
    siz = size(A);
    n = siz(1);
    A = reshape(A,n,[]);
    J = reshape(J,n,[]);
    
    allnan = all(J,2);
    i1 = find(~allnan,1)-1;
    i2 = find(~allnan,1,'last')+1;
    
    A([1:i1 i2:n],:) = [];
    A = reshape(A,[i2-i1-1 siz(2:end)]);
    A = shiftdim(A,1);
    
    I([1:i1 i2:n],:) = false;
    I = shiftdim(I,1);
    J = shiftdim(J,1);
end
end

%% Example
function RunTheExample
    % -- Cellular vortical flow --
    [x,y] = meshgrid(linspace(0,1,64));
    Vx = cos(2*pi*x+pi/2).*cos(2*pi*y);
    Vy = sin(2*pi*x+pi/2).*sin(2*pi*y);
    % -- Corrupt the original flow --
    Vx = Vx + sqrt(0.1)*randn(size(Vx)); % adding Gaussian noise
    Vy = Vy + sqrt(0.1)*randn(size(Vx)); % adding Gaussian noise
    I = randperm(numel(Vx));
    n = round(numel(Vx)/5);
    Vx(I(1:n)) = (rand(n,1)-0.5)*4; % adding outliers
    Vy(I(1:n)) = (rand(n,1)-0.5)*4; % adding outliers
    % -- Smooth the corrupt velocity field using PPPIV --
    wstate = warning('query','MATLAB:pppiv:LowSNR');
    warning('off','MATLAB:pppiv:LowSNR')
    [Vx2,Vy2] = pppiv(Vx,Vy);
    warning(wstate)
    % -- Display the results --
    figure, quiver(x,y,Vx,Vy,2.5), axis square
    title('Noisy velocity field')
    figure, quiver(x,y,Vx2,Vy2,2.5), axis square
    title('Post-processed velocity field')
end





