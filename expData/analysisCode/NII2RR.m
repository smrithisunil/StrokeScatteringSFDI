% Convert spectrum to reflectivity profile
% RAM = FRG raw file X 4

function RR = NII2RR(NII, intpDk, bavgfrm)

% choose parameter for lamda-k interpolation
if nargin < 2
	intpDk = -0.235;
end
if nargin < 3
	bavgfrm = 0;
end

% substract the reference signal, Subtract mean
[nk,nx,nf] = size(NII);
if bavgfrm == 1
    NII = NII - repmat(mean(NII(:,:),2),[1 nx nf]);
else
    for ifr=1:nf
        NII(:,:,ifr) = NII(:,:,ifr) - repmat(mean(NII(:,:,ifr),2),[1 nx]);
    end
end

% lamda-k interpolation
if intpDk ~= 0
    k = linspace(1-intpDk/2, 1+intpDk/2, nk);
    lam = 1./fliplr(k);
    for ifr=1:nf
        NII(:,:,ifr) = interp1(lam, NII(:,:,ifr), linspace(min(lam),max(lam),length(lam)), 'spline');
        if (mod(ifr,ceil(nf/5)) == 0)  
            disp(['... NIItoRR ' num2str(ifr) '/' num2str(nf) '  ' datestr(now,'HH:MM')]);  
        end
    end	
end

% test 45nm bandwidth light source
% hc = hamming(278); 
% h = [repmat(hc(1),(1024-278)/2,1); hc; repmat(hc(end),(1024-278)/2,1)];
% [nk,nx,ny]=size(NII);
% for ix=1:nx
%     for iy = 1:ny 
%         NII(:,ix,iy) = squeeze(NII(:,ix,iy)).*h;
%     end   
% end
    
% ifft
nz = round(nk/2);
RR = zeros(nz,nx,nf,'single');
for ifr=1:nf
    RRy = ifft(NII(:,:,ifr));
    RR(:,:,ifr) = RRy(1:nz,:);
end
