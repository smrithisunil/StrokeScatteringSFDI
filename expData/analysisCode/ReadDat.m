% load 3D data, Return nk_Nx_ny, applicable for non-stop Ascan
function [data]= ReadDat(filePath, dim, iseg) % [nk nx nframe/ny]
if nargin <3
    iseg=0;
end

nk = dim(1); nxRpt = dim(2); nx=dim(3); nyRpt = dim(4); ny = dim(5);
data = zeros(nk,nxRpt*nx*nyRpt,ny, 'single');
% read data   
fid=fopen(filePath,'r','l');
Start_iseg=(iseg-1)*nk*nxRpt*nx*nyRpt*ny*2;
for i = 1:ny
    fseek(fid,Start_iseg+(i-1)*nk*nxRpt*nx*nyRpt*2,'bof'); % due to the data type is int16, it therefore has to x2
    frame_data = fread(fid, nk*nxRpt*nx*nyRpt, 'int16');
    data(:,:,i) = reshape(frame_data, [nk nxRpt*nx*nyRpt]);
    if (mod(i,ceil(ny/2)) == 0)  
        disp(['... ReadDat ' num2str(i) '/' num2str(ny) '	' datestr(now,'HH:MM')]);  
    end    
end
fclose(fid);
   

   