function SS = gaussianss_dog(I,min_scale,max_scale,S,sigma0)

if(~isreal(I) || ndims(I) > 2)
    error('I must be a real two dimensional matrix') ;
end


scales=sigma0 * 2^(1/S) .^ (0:100);
flg = find(min_scale<=scales & scales<=max_scale);
flg = [flg(1)-2 flg(1)-1 flg flg(end)+1];
flg=flg(1<=flg);
SS.scales=scales(flg);

ftransform = fft2(I);
[ysize xsize] = size(ftransform);
[x y] = meshgrid((0 : xsize-1)-(xsize-1)/2, (0 : ysize-1)-(ysize-1)/2);
for i=1:length(SS.scales)
    s = SS.scales(i);
    SS.octave(:,:,i) = ifftshift(ifft2(fft2(exp(-(x .^ 2 + y .^ 2) / s ^ 2 /2) / pi / s ^ 2 /2) .* ftransform));
    SS.Lxx(:,:,i)= subdxx(SS.octave(:,:,i),1);
    SS.Lyy(:,:,i)= subdyy(SS.octave(:,:,i),1);
    SS.Lxy(:,:,i)= subdxy(SS.octave(:,:,i),1);
end
for i=1:length(SS.scales)-1
    SS.doG(:,:,i)= SS.scales(i)*(SS.octave(:,:,i+1) - SS.octave(:,:,i))/(SS.scales(i+1)-SS.scales(i));
end