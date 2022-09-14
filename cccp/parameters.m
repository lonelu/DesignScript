chains = 4;
chL = 28;
r0 = 7.36;
r1 = 2.26;
w0 = -2.45;
w1 = 102.857;
a = -12.01;
ph1 = [-9 -9 -9 -9];
cr = [0, 1, 0];
zoff = [0 0 0];
dph0 = [90 180 270];
%varargin.zoffaa = 1;
varargin.registerzoff = 1;
%varargin.apNNzoff = 1;
%varargin = {}

XYZ = generateCrickBB(chains, chL, r0, r1, w0, w1, a, ph1, cr, dph0, zoff, varargin);

%writematrix(XYZ,'XYZ.txt','Delimiter','tab')