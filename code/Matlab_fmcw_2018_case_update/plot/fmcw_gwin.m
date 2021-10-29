function h = fmcw_gwin(winfunh,varargin)

% fmcw_gblack(winfunh,varargin)
%
% args: winfunh funtion handle for window (e.g. @blackman)
% varargin: plot args for line
% 
% Plot a window centred at the mouse click (in dB scale)
% note: this is only valid for standard centre ferquency and bandwidth
%
% Craig Stewart
% 2014/6/5
%
% example use
%h(1) = fmcw_gwin(@blackman,'r');
%h(2) = fmcw_gwin(@blackmanharris,'g');
%h(3) = fmcw_gwin(@nuttallwin,'b');
%legend(h,{'blackman','blackmanharris','nuttallwin'})

if nargin == 0
    winfunh = @blackman;
end

% Plot parameters
p = 10; % pad factor (interpolate)
n = 99; % number of frequency bins to plot (must be odd)
N = n*p; % total points

% Radar parameters
er = 3.1;
c = 3e8; % velocity
B = 2e8; % bandwidth
dr = c/(2*p*B*sqrt(er)); % Brennan eq 14 rearranged 

% Get centre point
[x0,y0] = ginput(1);
if isempty(x0)
    x0 = 0;
    y0 = 0;
end

% Make window (zero padding and circshifting by half win length to keep phase flat)
w  = window(winfunh,n);
wpad = [w; zeros(n*(p-1),1)];
xn = 0.5*(n-1);
wpads = circshift(wpad,[-xn 0]); % phase centre
a = fft(wpads);
a = fftshift(a);
a = a(2:end);

% Recentre at mouse click
x = x0 + dr*[-N/2+1:N/2-1];
y = y0 + 20*log10(abs(a))- 20*log10(max(abs(a)));

% Plot
% Amplitude
hold on,
h = plot(x,y,varargin{:});

% % Phase (if a phase axis is present)
% pax = findobj('tag','pax');
% if ~isempty(pax);
%     axes(pax)
%     hold on
%     h(2) = plot(x,angle(a),varargin{:});
% end
