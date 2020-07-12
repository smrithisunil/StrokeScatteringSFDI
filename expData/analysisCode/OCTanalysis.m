%% OCT analysis

%% read raw data
datapath  = 'I:\Scat7_SS55\OCT\1hr\';
filename='RAW-2048-0001-00500-002-1000-'; % filename without sequence
N_file_load=1; % total number of files to be loaded
int=0;
RR=0;
for file_i =3
    filePath=[datapath,filename,num2str(file_i),'.dat'];
    iseg=0;
    disp (['Now processing file:' num2str(file_i) '/' num2str(N_file_load) ]);
    dim= [2048 0001 0500 002 1000];
    [data]= ReadDat(filePath, dim, iseg);
    RR = NII2RR(data);
    [nz,nx,ny] = size(RR);
    RR1(:,:,:,1) = RR(:,1:nx/2,:);
    RR1(:,:,:,2) = RR(:,nx/2+1:end,:);
    clear RR
    I1=abs(RR1(:,:,:,1));
    I2=abs(RR1(:,:,:,2));
    clear RR1 
    I1=I1.^2;
    I2=I2.^2;
    I=(I1+I2)/2;
    I =log(I);
    int=I+int;
    clear I I1 I2
end
int=int./N_file_load;
int = permute(int,[2 3 1]);
int = imresize(int,0.5,'bicubic');
%% int profile in single plane
figure
img3= squeeze(int(125,:,1:400)); % z,y,x
h3=imshow(img3,[-6 0]); colormap gray  %default [-6 0], [0 1] for visible light
h3=gca;
h3.Visible='On';
%% slope map of signal in z
f=waitbar(0, [{'Please wait...'}, 'Now processing:' '0'], 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);
z=(170:320); % z-level for profiling, 150px for 5x, 50px for 10x
   z=z';
   Z_mm=linspace(0,3.3*length(z)*1e-3,length(z))';
K = zeros(250,500);

for x= 1:1:500  %dimension x of I 
    for y= 1:1:250  %dimension y of I
        temp2 =polyfit(Z_mm,squeeze((int(y,x,z))),1);
        K(y,x)= temp2(1);
    end
    waitbar(x/500,f, [{'Please wait...'}, 'Now processing:' num2str(x)], 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    if getappdata(f,'canceling')
        break
    end 
end

delete(f)
figure
imshow(K, 'DisplayRange', [-7 0]), colormap gray, colorbar % Default [-7 0] for 5x.
save([datapath,'\','OCTattenuation.mat'], 'K')
