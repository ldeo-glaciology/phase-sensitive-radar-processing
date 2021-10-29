function h = ershade(x,eh,el,patchcol,varargin)

% h = ershade(x,eh,el,patchcol,varargin)
%
% simple error shading, inputs must be vectors
%
% Craig Stewart
% 2014/11/16

if nargin<3
    el = -eh; % symetrical errors
end
if nargin<4
    patchcol = [0.9 0.9 0.9];
end
x = x(:)';
eh = eh(:)';
el = el(:)';
if length(eh) == 1
    eh = repmat(eh,size(x));
end
if length(el) == 1
    el = repmat(el,size(x));
end
nx = numel(x);
neh = numel(eh);
nel = numel(el);
if std([nx neh nel])
    error('inputs must be the same size')
end

% Plot
h = patch([x fliplr(x)],[eh fliplr(el)],patchcol,'edgeColor','none',varargin{:}); % 