clear;
I1=imreadbw('cameraman.tif') ;
figure;
imshow(I1);

I1=I1-min(I1(:)) ;
I1=I1/max(I1(:)) ;
[frames1] = doh_detector( I1) ;
% [frames1] = log_detector(I1) ;
hold on;
for ii=1:size(frames1,2)
    x0=frames1(1,ii);
    y0=frames1(2,ii);
    aa =frames1(8,ii);
    s1 = frames1(4,ii);
    s2 = frames1(5,ii);
    c1 = frames1(6,ii);
    d1 = frames1(7,ii);
    rmat=[cos(aa) -sin(aa); sin(aa) cos(aa)]';
    sigma1=transpose(rmat)*[s1^2 0; 0 s2^2]*rmat;
    sigma=sqrtm(sigma1);
    X = sigma*[cos(linspace(0,2*pi,30)) ; sin(linspace(0,2*pi,30)) ;] ;
    X(1,:) = X(1,:) + x0;
    X(2,:) = X(2,:) + y0;
    plot(X(1,:),X(2,:),'g');
end