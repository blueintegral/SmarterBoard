function oframes1=dog_detector_3D(I,varargin)

[M,N,T] = size(I) ;
sigma0 = 0.5;
presmooth = 0;
thresh =  0.025;
asp =  10;
iso=0;
S      =  3 ;
min_scale= 0;
max_scale      =  min(size(I))/6;
slow = 0;
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
        case 'slow'
            slow = varargin{k+1} ;
        otherwise
            error(['Unknown parameter ''' varargin{k} '''.']) ;
    end
end
% 
% if C > 1
%     error('I should be a grayscale image') ;
% end
% if presmooth~=0
%     I = imsmooth1(I,presmooth);
% end
gss = gaussianss_dog_3D(I,min_scale,max_scale,S,sigma0);
[M,N,T,S]=size(gss.doG);
idx = localmax_3D(gss.doG, 1e-7) ;
idx = [idx localmax_3D(-gss.doG,-1e-7)];
[i,j,k,s] = ind2sub( size( gss.doG ), idx ) ;
oframes = [j(:)';i(:)';k(:)';s(:)'] ;
oframes(1:4,:) = oframes(1:4,:) - 1;
oframes = refinemx_3D(oframes,gss.doG) ;
oframes(1:4,:) = oframes(1:4,:) + 1;
oframes = oframes(:,oframes(1,:)>1 & oframes(2,:)>1 & oframes(3,:)>1 & oframes(4,:)>1 &oframes(1,:)<M & oframes(2,:)<N & oframes(3,:)<T & oframes(4,:)<S);
cord0=floor(oframes(1:4,:));
cord1=ceil(oframes(1:4,:));
cordr=round(oframes(1:4,:));
x=oframes(1,:);
y=oframes(2,:);
z=oframes(3,:);
s=oframes(4,:);
x1=cord0(1,:);
y1=cord0(2,:);
z1=cord0(3,:);
s1=cord0(4,:);
x2=cord1(1,:);
y2=cord1(2,:);
z2=cord1(3,:);
s2=cord1(4,:);
xr=cordr(1,:);
yr=cordr(2,:);
zr=cordr(3,:);
sr=cordr(4,:);
ind1=sub2ind(size(gss.doG),y1,x1,z1,s1);
ind2=sub2ind(size(gss.doG),y1,x1,z2,s1);
ind3=sub2ind(size(gss.doG),y2,x1,z1,s1);
ind4=sub2ind(size(gss.doG),y2,x1,z2,s1);
ind5=sub2ind(size(gss.doG),y1,x2,z1,s1);
ind6=sub2ind(size(gss.doG),y1,x2,z2,s1);
ind7=sub2ind(size(gss.doG),y2,x2,z1,s1);
ind8=sub2ind(size(gss.doG),y2,x2,z2,s1);
ind9=sub2ind(size(gss.doG),y1,x1,z1,s2);
ind10=sub2ind(size(gss.doG),y1,x1,z2,s2);
ind11=sub2ind(size(gss.doG),y2,x1,z1,s2);
ind12=sub2ind(size(gss.doG),y2,x1,z2,s2);
ind13=sub2ind(size(gss.doG),y1,x2,z1,s2);
ind14=sub2ind(size(gss.doG),y1,x2,z2,s2);
ind15=sub2ind(size(gss.doG),y2,x2,z1,s2);
ind16=sub2ind(size(gss.doG),y2,x2,z2,s2);

indr=sub2ind(size(gss.doG),yr,xr,zr,sr); %20130606
x20=x2-x;
y20=y2-y;
z20=z2-z;
s20=s2-s;
x01=x-x1;
y01=y-y1;
z01=z-z1;
s01=s-s1;
x21=x2-x1;
y21=y2-y1;
z21=z2-z1;
s21=s2-s1;

m1= x20 .* y20 .* z20.*s20 ./ x21 ./ y21 ./ z21./s21 ;
m2= x20 .* y20 .* z01.*s20 ./ x21 ./ y21 ./ z21./s21 ;
m3= x20 .* y01 .* z20.*s20 ./ x21 ./ y21 ./ z21./s21 ;
m4= x20 .* y01 .* z01.*s20 ./ x21 ./ y21 ./ z21./s21 ;
m5= x01 .* y20 .* z20.*s20 ./ x21 ./ y21 ./ z21./s21 ;
m6= x01 .* y20 .* z01.*s20 ./ x21 ./ y21 ./ z21./s21 ;
m7= x01 .* y01 .* z20.*s20 ./ x21 ./ y21 ./ z21./s21 ;
m8= x01 .* y01 .* z01.*s20 ./ x21 ./ y21 ./ z21./s21;
m9= x20 .* y20 .* z20.*s01 ./ x21 ./ y21 ./ z21./s21 ;
m10= x20 .* y20 .* z01.*s01 ./ x21 ./ y21 ./ z21./s21 ;
m11= x20 .* y01 .* z20.*s01 ./ x21 ./ y21 ./ z21./s21 ;
m12= x20 .* y01 .* z01.*s01 ./ x21 ./ y21 ./ z21./s21 ;
m13= x01 .* y20 .* z20.*s01 ./ x21 ./ y21 ./ z21./s21 ;
m14= x01 .* y20 .* z01.*s01 ./ x21 ./ y21 ./ z21./s21 ;
m15= x01 .* y01 .* z20.*s01 ./ x21 ./ y21 ./ z21./s21 ;
m16= x01 .* y01 .* z01.*s01 ./ x21 ./ y21 ./ z21./s21;

Lxx=gss.Lxx(ind1).*m1 + gss.Lxx(ind2) .* m2 + gss.Lxx(ind3) .* m3 + gss.Lxx(ind4) .* m4 + ...
    gss.Lxx(ind5) .* m5 + gss.Lxx(ind6) .* m6 + gss.Lxx(ind7) .* m7 + gss.Lxx(ind8) .* m8 + ...
    gss.Lxx(ind9) .* m9 + gss.Lxx(ind10) .* m10 + gss.Lxx(ind11) .* m11 + gss.Lxx(ind12) .* m12 + ...
    gss.Lxx(ind13) .* m13 + gss.Lxx(ind14) .* m14 + gss.Lxx(ind15) .* m15 + gss.Lxx(ind16) .* m16;
nanidx =isnan(Lxx); 
Lxx(nanidx)=gss.Lxx(indr(nanidx));

Lyy=gss.Lyy(ind1).*m1 + gss.Lyy(ind2) .* m2 + gss.Lyy(ind3) .* m3 + gss.Lyy(ind4) .* m4 + ...
    gss.Lyy(ind5) .* m5 + gss.Lyy(ind6) .* m6 + gss.Lyy(ind7) .* m7 + gss.Lyy(ind8) .* m8 + ...
    gss.Lyy(ind9) .* m9 + gss.Lyy(ind10) .* m10 + gss.Lyy(ind11) .* m11 + gss.Lyy(ind12) .* m12 + ...
    gss.Lyy(ind13) .* m13 + gss.Lyy(ind14) .* m14 + gss.Lyy(ind15) .* m15 + gss.Lyy(ind16) .* m16;
Lyy(nanidx)=gss.Lyy(indr(nanidx));

Lzz=gss.Lzz(ind1).*m1 + gss.Lzz(ind2) .* m2 + gss.Lzz(ind3) .* m3 + gss.Lzz(ind4) .* m4 + ...
    gss.Lzz(ind5) .* m5 + gss.Lzz(ind6) .* m6 + gss.Lzz(ind7) .* m7 + gss.Lzz(ind8) .* m8 + ...
    gss.Lzz(ind9) .* m9 + gss.Lzz(ind10) .* m10 + gss.Lzz(ind11) .* m11 + gss.Lzz(ind12) .* m12 + ...
    gss.Lzz(ind13) .* m13 + gss.Lzz(ind14) .* m14 + gss.Lzz(ind15) .* m15 + gss.Lzz(ind16) .* m16;
Lzz(nanidx)=gss.Lzz(indr(nanidx));

Lxy=gss.Lxy(ind1).*m1 + gss.Lxy(ind2) .* m2 + gss.Lxy(ind3) .* m3 + gss.Lxy(ind4) .* m4 + ...
    gss.Lxy(ind5) .* m5 + gss.Lxy(ind6) .* m6 + gss.Lxy(ind7) .* m7 + gss.Lxy(ind8) .* m8 + ...
    gss.Lxy(ind9) .* m9 + gss.Lxy(ind10) .* m10 + gss.Lxy(ind11) .* m11 + gss.Lxy(ind12) .* m12 + ...
    gss.Lxy(ind13) .* m13 + gss.Lxy(ind14) .* m14 + gss.Lxy(ind15) .* m15 + gss.Lxy(ind16) .* m16;
Lxy(nanidx)=gss.Lxy(indr(nanidx));

Lxz=gss.Lxz(ind1).*m1 + gss.Lxz(ind2) .* m2 + gss.Lxz(ind3) .* m3 + gss.Lxz(ind4) .* m4 + ...
    gss.Lxz(ind5) .* m5 + gss.Lxz(ind6) .* m6 + gss.Lxz(ind7) .* m7 + gss.Lxz(ind8) .* m8 + ...
    gss.Lxz(ind9) .* m9 + gss.Lxz(ind10) .* m10 + gss.Lxz(ind11) .* m11 + gss.Lxz(ind12) .* m12 + ...
    gss.Lxz(ind13) .* m13 + gss.Lxz(ind14) .* m14 + gss.Lxz(ind15) .* m15 + gss.Lxz(ind16) .* m16;
Lxz(nanidx)=gss.Lxz(indr(nanidx));

Lyz=gss.Lyz(ind1).*m1 + gss.Lyz(ind2) .* m2 + gss.Lyz(ind3) .* m3 + gss.Lyz(ind4) .* m4 + ...
    gss.Lyz(ind5) .* m5 + gss.Lyz(ind6) .* m6 + gss.Lyz(ind7) .* m7 + gss.Lyz(ind8) .* m8 + ...
    gss.Lyz(ind9) .* m9 + gss.Lyz(ind10) .* m10 + gss.Lyz(ind11) .* m11 + gss.Lyz(ind12) .* m12 + ...
    gss.Lyz(ind13) .* m13 + gss.Lyz(ind14) .* m14 + gss.Lyz(ind15) .* m15 + gss.Lyz(ind16) .* m16;
Lyz(nanidx)=gss.Lyz(indr(nanidx));

doG=gss.doG(ind1).*m1 + gss.doG(ind2) .* m2 + gss.doG(ind3) .* m3 + gss.doG(ind4) .* m4 + ...
    gss.doG(ind5) .* m5 + gss.doG(ind6) .* m6 + gss.doG(ind7) .* m7 + gss.doG(ind8) .* m8 + ...
    gss.doG(ind9) .* m9 + gss.doG(ind10) .* m10 + gss.doG(ind11) .* m11 + gss.doG(ind12) .* m12 + ...
    gss.doG(ind13) .* m13 + gss.doG(ind14) .* m14 + gss.doG(ind15) .* m15 + gss.doG(ind16) .* m16;
doG(nanidx)=gss.doG(indr(nanidx));

L=gss.L(ind1).*m1 + gss.L(ind2) .* m2 + gss.L(ind3) .* m3 + gss.L(ind4) .* m4 + ...
    gss.L(ind5) .* m5 + gss.L(ind6) .* m6 + gss.L(ind7) .* m7 + gss.L(ind8) .* m8 + ...
    gss.L(ind9) .* m9 + gss.L(ind10) .* m10 + gss.L(ind11) .* m11 + gss.L(ind12) .* m12 + ...
    gss.L(ind13) .* m13 + gss.L(ind14) .* m14 + gss.L(ind15) .* m15 + gss.L(ind16) .* m16;
L(nanidx)=gss.L(indr(nanidx));

% %noise suppression
% tt=[gss.doG(ind1);gss.doG(ind2);gss.doG(ind3);gss.doG(ind4);gss.doG(ind5);gss.doG(ind6);gss.doG(ind7);gss.doG(ind8)];
% f0 = sum(tt>0)>=5;
% tt(tt<0)=nan;
% mt=min(tt);
% Mt=max(tt);
% f1= f0 & 0.88<mt./Mt;
% 
oframes1=[];
for ii=1:length(Lxx)
    hess = [Lxx(ii) Lxy(ii) Lxz(ii);Lxy(ii) Lyy(ii) Lyz(ii);Lxz(ii) Lyz(ii) Lzz(ii)];
    [V,D]=eig(hess);
    D=diag(D);
    if sum(D<0)~=3 && sum(D>0)~=3
        continue;
    end
    [B,IX]=sort(abs(D),'descend');
    D=D(IX);
    V= V(:,IX);
    R1=D(1)/D(2);
    R2=D(1)/D(3);
    K=R1 * (R2 ^ 2 + R1 ^ 2 * (3 + 2 * R2 + 3 * R2 ^ 2)) / (2 * R1 * R2 + 3 * R2 ^ 2 + R1 ^ 2 * (3 + R2 ^ 2));
    G=R2 * (3 * R2 ^ 2 + 2 * R1 * R2 ^ 2 + R1 ^ 2 * (1 + 3 * R2 ^ 2)) / R1 / (R2 ^ 2 + R1 ^ 2 * (3 + 2 * R2 + 3 * R2 ^ 2));
    H=((2 * R1 * R2 + 3 * R2 ^ 2 + R1 ^ 2 * (3 + R2 ^ 2)) / R1 / R2 / (R1 + R2 + R1 * R2)) / 0.2e1;
    if K<0 ||G<0||H<0
        continue;
    end
    c =  -doG(ii) * G ^ (-0.1e1 / 0.2e1) / K * (H / (1 + H) / (1 + H * K) / (1 + G * H * K)) ^ (-0.3e1 / 0.2e1) / (3 + H ^ 2 * K * (1 + G + G * K) + 2 * H * (1 + K + G * K));
    if abs(c)<thresh
        continue;
    end
    d =-c * G * H ^ (0.3e1 / 0.2e1) * K * (G * (1 + H) * (1 + H * K) * (1 + G * H * K)) ^ (-0.1e1 / 0.2e1) + L(ii);
    oframes1=[oframes1 [oframes(:,ii);zeros(size(K));H;K;G;c;d;V(:,1)]];
end
oframes = [];
if isempty(oframes1)
    return;
end
oframes1(5,:)=interp1(1:length(gss.scales),gss.scales,oframes1(4,:),'spline');
oframes1(6,:) = oframes1(5,:).*sqrt(oframes1(6,:));
oframes1(7,:) = oframes1(6,:).*sqrt(oframes1(7,:));
oframes1(8,:) = oframes1(7,:).*sqrt(oframes1(8,:));
end