clear;
NOISY = 0;
%download a external routine before setting EXTERNAL to 1
EXTERNAL = 0;
ISO=0;

%-pi/2 pi/2
angle = -pi/2 + pi.*rand(1000,1);
%2-15
radius = 2  + 13.*rand(1000,1);
baseline=rand(1000,1);
%1,9
aspect_ratio = 1 + 8.*rand(1000,1);
x = -40 + 80.*rand(1000,1);
y = -40 + 80.*rand(1000,1);
all = [];
all_asp=[];
WINSIZ=127;
k = 1;

for j=1:1000
    fprintf('%d.',j);
    d = baseline(j);
    c = -d + rand();
    ang = angle(j);
    norm_r = radius(j);
    asp = aspect_ratio(j);
    I1=gauss2(ang,norm_r/sqrt(asp),norm_r*sqrt(asp),WINSIZ,x(j),y(j));
    I1=I1-min(I1(:)) ;
    I1=I1/max(I1(:)) ;
    I1=I1*c + d ;
    if NOISY
        I1= imnoise(I1,'gaussian');
    end
    if mod(j,2)==1
        I1(end,:) = [];
        I1(:,end) = [];
    end
    figure(1);
    clf;
    subplot(1,2,1);
    imshow(I1);
    hold on;
    s1 = norm_r/sqrt(asp);
    s2 = norm_r*sqrt(asp);
    rmat=[cos(ang) sin(ang); -sin(ang) cos(ang)];
    sigma1=transpose(rmat)*[s1^2 0; 0 s2^2]*rmat;
    sigma=sqrtm(sigma1);
    X = sigma*[cos(linspace(0,2*pi,30)) ; sin(linspace(0,2*pi,30)) ;] ;
    X(1,:) = X(1,:) + x(j)+128;
    X(2,:) = X(2,:) + y(j)+128;
    hold on;
    line(X(1,:),X(2,:),'Color','g','LineWidth',1);
    title('expected');
    if EXTERNAL
        [frames1] = doh_detector( I1,'file','I1.dohaff','iso',ISO);
    else
        [frames1] = doh_detector( I1) ;
        %         [frames1] = log_detector( I1) ;
        %         [frames1] = dog_detector( I1) ;
        %         [frames1] = doh_detector_II( I1) ;
        %         [frames1] = log_detector_II( I1) ;
        
    end
    if isempty(frames1)
        continue;
    end
    
    subplot(1,2,2);
    if EXTERNAL
        imwrite(I1,'I1.ppm');  display_features('I1.dohaff','I1.ppm',1,1);
    else
        imshow(I1);
        hold on;
        for ii=1:size(frames1,2)
            x0=frames1(1,ii);
            y0=frames1(2,ii);
            if ~NOISY
                fprintf(' %f\t%f\t%f\t%f\t%f\t%f\t%f\n',x0-128-x(j), y0-128-y(j),100*(frames1(6,ii)-c)./c, 100*(frames1(7,ii)-d)./d, 100*(frames1(4,ii)-norm_r/sqrt(asp))./(norm_r/sqrt(asp)),100*(frames1(5,ii)-norm_r*sqrt(asp))./(norm_r*sqrt(asp)),100*(frames1(8,ii)-ang)./ang);
            end
            if ISO
                s1 = k*frames1(9,:);
                s2 = k*frames1(9,:);
                aa = 0;
            else
                aa = frames1(8,ii);
                s1 = k*frames1(4,ii);
                s2 = k*frames1(5,ii);
            end
            rmat=[cos(aa) sin(aa); -sin(aa) cos(aa)];
            sigma1=transpose(rmat)*[s1^2 0; 0 s2^2]*rmat;
            sigma=sqrtm(sigma1);
            X = sigma*[cos(linspace(0,2*pi,30)) ; sin(linspace(0,2*pi,30)) ;] ;
            X(1,:) = X(1,:) + x0;
            X(2,:) = X(2,:) + y0;
            hold on;
            line(X(1,:),X(2,:),'Color','g','LineWidth',1);
        end
    end
    drawnow;
    title('detected');
    pause(.1);
end