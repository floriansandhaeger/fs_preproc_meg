function source = fs_beamform(cfg,C,leadfield)

% USE DUMMY COVARIANCE MATRIX
if isempty(C)	
	C_all = zeros(length(leadfield.cfg.channel));
	for sidx = 1:length(leadfield.leadfield),
		if ~isnan(leadfield.leadfield{sidx}),
			C_act = leadfield.leadfield{sidx}*leadfield.leadfield{sidx}';
			C_all = C_all + C_act;
		end
	end
	C = C_all;
end

% regularization factor
if ~isfield(cfg,'a'), cfg.a = 0.05; end; a = cfg.a;

if ~isfield(cfg,'pow'), cfg.pow = 0; end;

if ~isfield(cfg,'combine'), cfg.combine = 0; end % 0: don't combine. 1: max variance

if islogical(leadfield.inside),
	leadfield.inside = find(leadfield.inside);
end

% compute regularization factor
a = a * trace(C)/size(C,1);
%regularize and invert CSD
Cinv = pinv(real(C) + a * eye(size(C)));
% imaginary part doesn't contain volume conduction, and for beamforming we're interested in that!

%compute noise level estimate
noise = svd(C); noise = noise(end); noise = max(noise,a);
			
clear source;

source.pow = repmat((NaN),size(leadfield.leadfield,1),size(leadfield.leadfield,2));
source.noise = repmat((NaN),size(leadfield.leadfield,1),size(leadfield.leadfield,2));
source.filt_raw = nan(3,size(C,1),length(leadfield.inside));
Creal = real(C);

for s=1:length(leadfield.inside),
	
	L = leadfield.leadfield{leadfield.inside(s)};
	filter = pinv(L' * Cinv * L) * L' * Cinv;
	CSD{s} = filter * Creal * ctranspose(filter);
	source.filt_raw(:,:,s) = filter;
	
	if cfg.combine,
		[so, ori] = lambda1(CSD{s});
		source.ori(s,:) = ori;
		filt_comb(:,s) = ((squeeze(source.filt_raw(:,:,s)))'*(source.ori(s,:))')';
	end
	
	if cfg.pow,
		% A)
		[so, ori] = lambda1(CSD{s});
		source.pow(leadfield.inside(s)) = so;
		source.noise(leadfield.inside(s)) = noise * lambda1(filter * ctranspose(filter));
		% B)
		%source.pow(leadfield.inside(s)) = real(trace(CSD{s}));
		%source.noise(leadfield.inside(s)) = noise * real(trace(source.filter{s} * ctranspose(source.filter{s})));
		
		%source.filt_comb(s,:) = (source.filter{s}'*ori)';                          % maps on dominant direction, will be a voxels * sensors matrix
		
		source.ori(s,:) = ori;
		
		filt_comb(:,s) = ((squeeze(source.filt_raw(:,:,s)))'*(source.ori(s,:))')';
	end
end

if cfg.combine,
	source.filt_comb = filt_comb;
end

if cfg.pow,
	noise = source.noise(~isnan(source.noise));
	pow = source.pow(~isnan(source.pow));
	source.nai = pow./noise;
	source.pow = pow;
	source.noise = noise;
	source.filt_comb = filt_comb;
end

source.C       = C;
source.pos     = leadfield.pos;
source.inside  = leadfield.inside;
source.outside = leadfield.outside;
source.unit    = leadfield.unit;
if isfield(leadfield,'dim'), source.dim   = leadfield.dim; end
%source.vol      = leadfieldfile;
%source.leadfield = leadfieldfile;



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute the pseudo inverse. This is the same as the
% standard Matlab function, except that the default tolerance is twice as
% high.
%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 8127 $  $Date: 2009/06/17 13:40:37 $
%   default tolerance increased by factor 2 (Robert Oostenveld, 7 Feb 2004)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = pinv(A,varargin)
[m,n] = size(A);
if n > m
    X = pinv(A',varargin{:})';
else
    [U,S,V] = svd(A,0);
    if m > 1, s = diag(S);
    elseif m == 1, s = S(1);
    else s = 0;
    end
    if nargin == 2
        tol = varargin{1};
    else
        tol = 10 * max(m,n) * max(s) * eps;
    end
    r = sum(s > tol);
    if (r == 0)
        X = zeros(size(A'),class(A));
    else
        s = diag(ones(r,1)./s(1:r));
        X = V(:,1:r)*s*U(:,1:r)';
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to obtain the largest singular value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s, ori] = lambda1(x)
% determine the largest singular value, which corresponds to the power along the dominant direction
[u, s, v] = svd(x);
s   = s(1);
ori = u(:,1);
end