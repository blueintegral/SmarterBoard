function SS = gaussianss_doh_II(I,min_scale,max_scale,S,sigma0)
if(~isreal(I) || ndims(I) > 2)
    error('I must be a real two dimensional matrix') ;
end
scales=sigma0 * 2^(1/S) .^ (0:10);
scales(scales<1) = [];
scales = round(scales);
if scales(1) == scales(2)
    scales(1)=[];
end
SS.scales=scales(min_scale<=scales & scales<=(max_scale));
II = cumsum(cumsum(I,2),1);
for i=1:length(SS.scales)
    s = SS.scales(i);
    s4 = SS.scales(i)^4;
    SS.L(:,:,i)  = sub4(II,s);
    SS.Lxx(:,:,i)= subdxx(SS.L(:,:,i),s);
    SS.Lyy(:,:,i)= subdyy(SS.L(:,:,i),s);
    SS.Lxy(:,:,i)= subdxy(SS.L(:,:,i),s);
    SS.doh(:,:,i)= ( SS.Lxx(:,:,i).*SS.Lyy(:,:,i)-SS.Lxy(:,:,i).*SS.Lxy(:,:,i))*s4;
end