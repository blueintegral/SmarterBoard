function [oframes]=log_detector_II(I,varargin)
% The lines of ¡±oframes¡± represents:
% 1 x
% 2 y
% 3 interpolated index of detection scale
% 4 width
% 5 length
% 6 contrast
% 7 offset
% 8 angle
% 9 detection scale
[M,N,C] = size(I) ;
S      =  2 ;
min_scale= 1;
max_scale = round(min([M,N])/7);
sigma0 = 1;                  % smooth lev. -1 at 1.6
presmooth = 0.75;
thresh =  0.025;
asp =  5;
iso=0;
for k=1:2:length(varargin)
    switch lower(varargin{k})
        case 'asp'
            asp = varargin{k+1} ;
            
        case 'min_scale'
            min_scale = varargin{k+1} ;
        case 'max_scale'
            max_scale = varargin{k+1} ;
            
        case 'numlevels'
            S = varargin{k+1} ;
            
        case 'sigma0'
            sigma0 = varargin{k+1} ;
            
        case 'presmooth'
            presmooth = varargin{k+1} ;
            
        case 'threshold'
            thresh = varargin{k+1} ;
        case 'file'
            file = varargin{k+1} ;
        case 'iso'
            iso = varargin{k+1} ;
        otherwise
            error(['Unknown parameter ''' varargin{k} '''.']) ;
    end
end
if mod(M,2) == 0
    I(end,:,:) = [];
end
if  mod(N,2) == 0
    I(:,end,:) = [];
end
if C > 1
    error('I should be a grayscale image') ;
end
if presmooth~=0
    I = imsmooth1(I,presmooth);
end
gss = gaussianss_log_II(I,min_scale,max_scale,S,sigma0) ;
%c=0.01,sqrt(K)=50
idx = siftlocalmax(  gss.LoG, 0.00385208  ) ;
idx = [idx , siftlocalmax( - gss.LoG, 0.00385208)] ;

[i,j,s] = ind2sub( size( gss.LoG ), idx ) ;
oframes = [j(:)';i(:)';s(:)'] ;
oframes(1:3,:) = oframes(1:3,:) - 1;
oframes = dohrefinemx(oframes,gss.LoG) ;
oframes(1:3,:) = oframes(1:3,:) + 1;
cord0=floor(oframes(1:3,:));
cord1=ceil(oframes(1:3,:));
cordr=round(oframes(1:3,:));
x=oframes(1,:);
y=oframes(2,:);
z=oframes(3,:);
x1=cord0(1,:);
y1=cord0(2,:);
z1=cord0(3,:);
x2=cord1(1,:);
y2=cord1(2,:);
z2=cord1(3,:);
xr=cordr(1,:);
yr=cordr(2,:);
zr=cordr(3,:);
ind1=sub2ind(size(gss.LoG),y1,x1,z1);
ind2=sub2ind(size(gss.LoG),y1,x1,z2);
ind3=sub2ind(size(gss.LoG),y2,x1,z1);
ind4=sub2ind(size(gss.LoG),y2,x1,z2);
ind5=sub2ind(size(gss.LoG),y1,x2,z1);
ind6=sub2ind(size(gss.LoG),y1,x2,z2);
ind7=sub2ind(size(gss.LoG),y2,x2,z1);
ind8=sub2ind(size(gss.LoG),y2,x2,z2);
indr=sub2ind(size(gss.LoG),yr,xr,zr); %20130606
x20=x2-x;
y20=y2-y;
z20=z2-z;
x01=x-x1;
y01=y-y1;
z01=z-z1;
x21=x2-x1;
y21=y2-y1;
z21=z2-z1;
Lxx=gss.Lxx(ind1) .* x20 .* y20 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxx(ind2) .* x20 .* y20 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxx(ind3) .* x20 .* y01 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxx(ind4) .* x20 .* y01 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxx(ind5) .* x01 .* y20 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxx(ind6) .* x01 .* y20 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxx(ind7) .* x01 .* y01 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxx(ind8) .* x01 .* y01 .* z01 ./ x21 ./ y21 ./ z21;
nanidx =isnan(Lxx); %20130606
Lxx(nanidx)=gss.Lxx(indr(nanidx));

Lxy=gss.Lxy(ind1) .* x20 .* y20 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxy(ind2) .* x20 .* y20 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxy(ind3) .* x20 .* y01 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxy(ind4) .* x20 .* y01 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxy(ind5) .* x01 .* y20 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxy(ind6) .* x01 .* y20 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxy(ind7) .* x01 .* y01 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.Lxy(ind8) .* x01 .* y01 .* z01 ./ x21 ./ y21 ./ z21;
Lxy(nanidx)=gss.Lxy(indr(nanidx));
Lyy=gss.Lyy(ind1) .* x20 .* y20 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.Lyy(ind2) .* x20 .* y20 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.Lyy(ind3) .* x20 .* y01 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.Lyy(ind4) .* x20 .* y01 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.Lyy(ind5) .* x01 .* y20 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.Lyy(ind6) .* x01 .* y20 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.Lyy(ind7) .* x01 .* y01 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.Lyy(ind8) .* x01 .* y01 .* z01 ./ x21 ./ y21 ./ z21;
Lyy(nanidx)=gss.Lyy(indr(nanidx));
L=gss.L(ind1) .* x20 .* y20 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.L(ind2) .* x20 .* y20 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.L(ind3) .* x20 .* y01 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.L(ind4) .* x20 .* y01 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.L(ind5) .* x01 .* y20 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.L(ind6) .* x01 .* y20 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.L(ind7) .* x01 .* y01 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.L(ind8) .* x01 .* y01 .* z01 ./ x21 ./ y21 ./ z21;
L(nanidx)=gss.L(indr(nanidx));
LoG=gss.LoG(ind1) .* x20 .* y20 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.LoG(ind2) .* x20 .* y20 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.LoG(ind3) .* x20 .* y01 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.LoG(ind4) .* x20 .* y01 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.LoG(ind5) .* x01 .* y20 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.LoG(ind6) .* x01 .* y20 .* z01 ./ x21 ./ y21 ./ z21 + ...
    gss.LoG(ind7) .* x01 .* y01 .* z20 ./ x21 ./ y21 ./ z21 + ...
    gss.LoG(ind8) .* x01 .* y01 .* z01 ./ x21 ./ y21 ./ z21;
LoG(nanidx)=gss.LoG(indr(nanidx));
tt=[gss.LoG(ind1);gss.LoG(ind2);gss.LoG(ind3);gss.LoG(ind4);gss.LoG(ind5);gss.LoG(ind6);gss.LoG(ind7);gss.LoG(ind8)];
mt=min(tt);
Mt=max(tt);
f1= 0.88<mt./Mt & mt./Mt<1/0.88;
Lcc= sqrt((4 * Lxy .^ 2 + Lxx .^ 2 - 2 .* Lxx .* Lyy + Lyy .^ 2));
e1=Lxx / 0.2e1 - Lcc / 0.2e1 + Lyy / 0.2e1;
e2= Lxx / 0.2e1 + Lcc / 0.2e1 + Lyy / 0.2e1;
r=e1./e2;
ind1 = r>0 & f1;
oframes=oframes(:,ind1);
if isempty(oframes)
    return;
end
Lxx = Lxx(ind1);
Lyy = Lyy(ind1);
Lxy = Lxy(ind1);
Lcc = Lcc(ind1);
L = L(ind1);
LoG= LoG(ind1);
r   = r(ind1);
flg = r<1;
r(flg) = 1./r(flg);
t1= Lxy * 0.2e1./(Lxx - Lcc - Lyy);
t2= Lxy * 0.2e1./(Lxx + Lcc - Lyy);
t1(flg) = t2(flg);
t1 = atan(t1);
K = (r + 3 * r .^ 3) ./ (3 + r .^ 2);
H = ((3 + r .^ 2) ./ r ./ (1 + r)) / 0.2e1;
c =   -(LoG.* (3 + 2 * r + 3 * r .^ 2) .^ 2 ./ (1 + r) .^ 2 .* (3 + 10 * r .^ 2 + 3 * r .^ 4) .^ (-0.1e1 / 0.2e1)) / 0.2e1/2;
d =  (-c .* sqrt((3 + 10 * r .^ 2 + 3 * r .^ 4)) + ((3 + 2 * r + 3 * r .^ 2) .* L)) ./ (3 + 2 * r + 3 * r .^ 2);
oframes = [oframes;sqrt(H);sqrt(K);c;d;t1;zeros(size(c))];
flg=oframes(5,:)<=asp & abs(oframes(6,:))>thresh;
oframes=oframes(:,flg);
if isempty(oframes)
    return
end
oframes(9,:)=interp1(1:length(gss.scales),gss.scales,oframes(3,:),'spline');
oframes(4,:) = oframes(4,:).*oframes(9,:);
oframes(5,:) = oframes(5,:).*oframes(4,:);
% resize proportionally as the external methods
% k = 3.5;
k = 1;
if iso
    a = k*oframes(9,:);
    b = k*oframes(9,:);
    theta = 0;
else
    a = k*oframes(4,:);
    b = k*oframes(5,:);
    theta = oframes(8,:);
end
u= (b .^ 2 .* cos(theta) .^ 2 + a .^ 2 .* sin(theta) .^ 2) ./ a .^ 2 ./ b .^ 2;
v= (b .^ 2 - a .^ 2) .* cos(theta) .* sin(theta) ./ a .^ 2 ./ b .^ 2;
w = (a .^ 2 .* cos(theta) .^ 2 + b .^ 2 .* sin(theta) .^ 2) ./ a .^ 2 ./ b .^ 2;
oframes = [oframes;u;v;w];
ybound=abs(u .* (u .* (-v .^ 2 + u .* w)) .^ (-0.1e1 / 0.2e1));
xbound = abs(sqrt(w) .* (-v .^ 2 + u .* w) .^ (-0.1e1 / 0.2e1));
sel = oframes(1,:)+xbound<N & oframes(1,:) - xbound>1 & oframes(2,:)+ ybound<M & oframes(2,:) - ybound>1 ;
oframes=oframes(:,sel) ;

if exist('file','var')
    f=fopen(file,'w');
    fprintf(f,'%d\n%d\n',0,size(oframes,2));
    for ii=1:size(oframes,2)
        fprintf(f,'%f %f %f %f %f\n',oframes(1,ii)-1,oframes(2,ii)-1,oframes(10,ii),oframes(11,ii),oframes(12,ii));
    end
    fclose(f);
end
end