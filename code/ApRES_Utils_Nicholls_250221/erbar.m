function [hb,hp] = erbar(x,y,eh,el,linecol,pointcol)

% [hb,hp] = erbar(x,y,eh,el,linecol,pointcol)
% simple errorbar, inputs must be vectors

if nargin<4
    el = -eh; % symetrical errors
end
if nargin<6
    pointcol = 'k';
    linecol = 'r';
end

x = reshape(x,1,numel(x));
y = reshape(y,1,numel(y));
eh = reshape(eh,1,numel(eh));
el = reshape(el,1,numel(el));
if length(eh) == 1
    eh = repmat(eh,size(x));
end
if length(el) == 1
    el = repmat(el,size(x));
end
% Plot
xb = [x; x];
yb = [y+eh; y+el];
hb = plot(xb,yb,'col',linecol,'linewidth',1.5); % Plot error bars
hold on
hp = plot(x,y,'.','col',pointcol); % Plot points