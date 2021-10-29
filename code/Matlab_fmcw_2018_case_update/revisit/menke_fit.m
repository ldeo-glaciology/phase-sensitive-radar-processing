function [M_est,M_error] = menke_fit(x,y,err)

% [M_est,M_error] = menke_fit(x,y,err)
% 
% Estimates the fit of the column vectors in 'x' to the data in 'y'
% and error of the the estimates based on the data error 'err'
% x and y with vector or scaler error estimates err using a weighted least
% squares approach ( see Geophysical Data Analysis: Discrete Inverse Theory
% revised edition by Wiliam Menke, Academic Press 1984 (1989) chapter 3)
%
% Note that this code can also be used for fitting any shape to the data
% and also for fitting planes etc.
%
% Note: Unless x contains a column of ones (or any other constant) then the
% fit will be forced through the origin
%
% Input: x = A Matrix of n column vectors each of which is a mode to fit to the data
%        (i.e. a column of ones gives the y intercept, a column of x values
%        gives the straight line slope...
%        y = vector of dependant variable
%        err = either (a) nothing (no errors),
%                    (b) a scalar representing the error on each y value
%                    (c) a vector of size(y) indicating the error on each y
%                    value
%
% Output: M_est = vector containing the estimated coefficients for each
%                column of x
%         M_error = the error associated with each M_est

% Program starts

% Test data
if nargin==0 % run test/demo
    figure
    for ii = 1:4
        % Linear
        slope = 2;
        offset = 2;
        er = ii*4; % random normal noise
        x = [0:100]';
        xrange = max(x)-min(x);
        yrange = slope*xrange;
        err = er*randn(size(x));
        y = offset + slope*x+err;
        X = [ones(size(x)) x]; % first column of ones fits a y intercept
        [M_est,M_error] = menke_fit(X,y,err);
        yest = M_est(1)+M_est(2)*x;
        
        % Plot
        ax(ii) = subplot(2,2,ii);
        title(['menke_fit demo: slope =' num2str(slope) ', err=' num2str(er)])
        datah = plot(x,y,'b.');
        hold on
        fith = plot(x,yest,'r');
        %legend([datah,fith],'data','fit')
        th = text(1,offset+20,{['fit: intercept =' num2str(M_est(1)) '+/-' num2str(num2str(M_error(1)))];...
            ['fit: slope     =' num2str(M_est(2)) '+/-' num2str(num2str(M_error(2)))]},'verticalAl','top');
        axis tight
        grid on
    end
    linkaxes(ax,'xy');
    return
end


%make sure the input data is a column vector
[m n] = size(x);
if (m < n), x = x'; end
[m n] = size(y);
if (m ==1), y = y'; end


%if no error vector is given then use an identity matrix
if (nargin == 2)
    W = diag(ones(1,length(x)));
else
    if length(err) == 1
        err = err*ones(size(y));
    end
    w = 1./(err.^2);
    W = diag(w);
end    

% Choose whether you want an intercept or to set it at zero
newcol = [];
% if all(std(x)>0);
%     answer = input('Do you want to force the fit through the origin? (y/n): ','s');
%     if strcmpi(answer,'n')
%         newcol = ones(size(x,1));
%     end
% end
G = [newcol x];

% Calculate the coefficients of each column of x
temp = inv(G'*W*G)* G' * W; % equation 3.37 page 54
M_est =  temp * y;          % equation 3.37 page 54

% Calculate the errors for the coefficients
if (nargin == 2)
    resid = y - G*M_est;
    err = cov(resid) * ((length(x)-1)/(max(size(G))-min(size(G)))) * inv(G' * G);
else
    w = err.^2;
    cov_y = diag(w);
    err = temp * cov_y * temp';
end
M_error = sqrt(diag(err));