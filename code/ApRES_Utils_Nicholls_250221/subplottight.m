function ah = subplottight(n,m,i,xmargin,ymargin)

% ah = subplottight(n,m,i)
%
% Subplot with no gaps between axes
% margin can be scalar (same in x and y) or a [xmargin ymargin] vector
% 
% probably from keith (due to lack of comments!) - rewritten to allow
% margins CLS 1 April 2014!

%[c,r] = ind2sub([m n], i);
%ah = subplot('Position', [(c-1)/m, 1-(r)/n, 1/m, 1/n]) % original ended
%here

[c,r] = ind2sub([m n], i);
if nargin == 3
    xmargin = 0.1;
    ymargin = 0.1;
end
if length(xmargin)==1
    xmargin = [xmargin xmargin];
end
if length(ymargin)==1
    ymargin = [ymargin ymargin];
end
w = (1-sum(xmargin))./m;
h = (1-sum(ymargin))./n;
x = xmargin(1) + w*(c-1);
y = 1-ymargin(2) - h*(r);
ah = subplot('Position', [x,y,w,h]);
%set(gca,'ActivePositionProperty','outerposition')