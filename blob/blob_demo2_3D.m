clear;
NOISY = 0;
%-pi/2 pi/2
angle = -pi/2 + pi.*rand(3,1000);
%2-15
r = 1  + 10.*rand(3,1000);
r = sort(r,1);
d=rand(1000,1);
c=-1+2*rand(1000,1);
ori = -10  + 20.*rand(3,1000);
W=30;
k = 1;

for j=1:1000
    fprintf('%d.',j);
    [I1,m1]=gauss3(r(:,j),W,ori(:,j),angle(:,j));
    %
    D=I1;
    D=uint8(D*255);
    figure(1);clf;
    Ds = smooth3(D);
    hiso = patch(isosurface(Ds,5),...
        'FaceColor',[1,.75,.65],...
        'EdgeColor','none');
    isonormals(Ds,hiso)
    hcap = patch(isocaps(D,5),...
        'FaceColor','interp',...
        'EdgeColor','none');
    view(35,30);
    L=1;
    H=2*W+1;
    axis([L H L H L H]);
    daspect([1,1,1])
    lightangle(45,30);
    set(gcf,'Renderer','zbuffer'); lighting phong
    set(hcap,'AmbientStrength',.6)
    set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50);
    %
    I1=I1-min(I1(:)) ;
    I1=I1/max(I1(:)) ;
    I1=I1*c(j) + d(j) ;
    [frames1] = doh_detector_3D( I1) ;
%     [frames1] = log_detector_3D( I1) ;
%     [frames1] = dog_detector_3D( I1) ;
    if isempty(frames1)
        continue;
    end
    for ii=1:size(frames1,2)
        x0=frames1(1,ii);
        y0=frames1(2,ii);
        z0=frames1(3,ii);
        if ~NOISY
            x=m1';y=frames1([11 12 13],ii);
            if y(1)<0
                y = -y;
            end
            angdiff=2 * atan(norm(x*norm(y) - norm(x)*y) / norm(x * norm(y) + norm(x) * y));
            fprintf(' %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',x0-W-1-ori(1,j), y0-W-1-ori(2,j),z0-W-1-ori(3,j),100*(frames1(9,ii)-c(j))./c(j), 100*(frames1(10,ii)-d(j))./d(j), 100*(frames1(6,ii)-r(1,j))./(r(1,j)),100*(frames1(7,ii)-r(2,j))./(r(2,j)),100*(frames1(8,ii)-r(3,j))./(r(3,j)),angdiff*180/pi);
        end
    end
end