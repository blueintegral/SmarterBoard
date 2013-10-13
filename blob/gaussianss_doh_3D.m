function SS = gaussianss_doh_3D(I,min_scale,max_scale,S,sigma0)
%
% if(~isreal(I) || ndims(I) > 2)
%     error('I must be a real two dimensional matrix') ;
% end

scales=sigma0 * 2^(1/S) .^ (0:100);
flg = find(min_scale<=scales & scales<=max_scale);
flg = [flg(1)-2 flg(1)-1 flg flg(end)+1];
flg=flg(1<=flg);
SS.scales=scales(flg);

ftransform = fftn(I);
[ysize xsize zsize] = size(ftransform);
[x y z] = meshgrid((0 : xsize-1)-(xsize-1)/2, (0 : ysize-1)-(ysize-1)/2,(0 : zsize-1)-(zsize-1)/2);
for i=1:length(SS.scales)
    s = SS.scales(i);
    s6 = s^6;
    gker=exp(-(x .^ 2 + y .^ 2 + z .^ 2) / s ^ 2 / 0.2e1) * sqrt(0.2e1) * pi ^ (-0.3e1 / 0.2e1) / s ^ 3 / 0.4e1;
    SS.L(:,:,:,i) = ifftshift(ifftn(fftn(gker) .* ftransform));
    Lxx=subdxx3(SS.L(:,:,:,i),1);
    Lyy=subdyy3(SS.L(:,:,:,i),1);
    Lzz=subdzz3(SS.L(:,:,:,i),1);
    Lxy=subdxy3(SS.L(:,:,:,i),1);
    Lxz=subdxz3(SS.L(:,:,:,i),1);
    Lyz=subdyz3(SS.L(:,:,:,i),1);
    SS.doH(:,:,:,i)= (-Lxz.^2.*Lyy+2*Lxy.*Lxz.*Lyz-Lxx.*Lyz.^2-Lxy.^2.*Lzz+Lxx.*Lyy.*Lzz)*s6;
    SS.Lxx(:,:,:,i) = Lxx;
    SS.Lyy(:,:,:,i) = Lyy;
    SS.Lzz(:,:,:,i) = Lzz;
    SS.Lxy(:,:,:,i) = Lxy;
    SS.Lxz(:,:,:,i) = Lxz;
    SS.Lyz(:,:,:,i) = Lyz;
end