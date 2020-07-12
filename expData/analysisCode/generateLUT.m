

function [reReflect, reReflect_dc] = generateLUT(fx, muaall, musall)

addpath('C:\BoasLab\Projects\SFDI\Code\MonteCarlo\Data')

Npoints=100; % number of points in the xy plane

x=(1:Npoints)*0.05;
kx=fx*2*pi;
Itemp1=cos(kx.*x)+1;
Itemp2=cos(kx.*x+(2/3*pi))+1;
Itemp3=cos(kx.*x+(4/3*pi))+1;

F = @(x) fftshift(fft2(ifftshift(x)));% define Fourier transform
Ft = @(x) fftshift(ifft2(ifftshift(x)));

load('IpointSource.mat')% for the phantom sample


Iillu2=repmat(Itemp1,[Npoints,1]);
Iillu2_2=repmat(Itemp2,[Npoints,1]);
Iillu2_3=repmat(Itemp3,[Npoints,1]);

If1p=Ft(F(Iillu2).*F(I));
If2p=Ft(F(Iillu2_2).*F(I));
If3p=Ft(F(Iillu2_3).*F(I));
MACp=(2^0.5)/3*sqrt((If1p-If2p).^2+(If2p-If3p).^2+(If3p-If1p).^2);
MDCp=1/3*(If1p+If2p+If3p);

for ii=1:length(muaall)
    for jj=1:length(musall)
        mus=musall(jj);
        mua=muaall(ii);
        filenm=['test_mus_',num2str(mus,'%.1f'),'_mua_',num2str(mua,'%.1f'),'.2ptout'];

        fid = fopen(filenm);
        Io = fread(fid,'float32');
        fclose(fid);
        I = reshape(Io,[Npoints Npoints]); 
        If1=Ft(F(Iillu2).*F(I));
        If2=Ft(F(Iillu2_2).*F(I));
        If3=Ft(F(Iillu2_3).*F(I));
        MAC=(2^0.5)/3*sqrt((If1-If2).^2+(If2-If3).^2+(If3-If1).^2);
        MDC=1/3*(If1+If2+If3);
        
        reReflect(ii,jj)=mean(mean(MAC./MACp));
        reReflect_dc(ii,jj)=mean(mean(MDC./MDCp));
        
    end
end

% imagesc(musall/10,muaall,reReflect);
% ylabel('\mu_a')
% xlabel('\mu_s''')
% colorbar
% title(['Relative Diffuse Reflectance f_x = ', num2str(fx), ' mm^{-1}'])

end