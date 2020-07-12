%% SFDI analysis

%%
trials = 5;
frequencies = 13;
phases = 3;
wavelengths = 1;
lambda = 530;
ref_exp = 50;
samp_exp = 200;
ref_exp = repmat(ref_exp, 1, phases*frequencies*trials);
samp_exp = repmat(samp_exp, 1, phases*frequencies*trials);

% Load reference sample 
pathname = uigetdir;
files = dir([pathname '\*.tif']);
for i = 1:trials*frequencies*phases*wavelengths
    data_ref(:,:,i) = double(imread([pathname '\' files(i).name]))./ref_exp(i); %.*power(i);   
end
data_ref = (data_ref);
data_ref = reshape(data_ref,[size(data_ref,1) size(data_ref,2) frequencies*phases*wavelengths trials]);
data_ref = mean(data_ref,4);

% Load the true sample
pathname = uigetdir;
files = dir([pathname '\*.tif']);
for i = 1:trials*frequencies*phases*wavelengths
    data(:,:,i) = double(imread([pathname '\' files(i).name]))./samp_exp(i); %.*power(i);
end
data = (data);
data = reshape(data,[size(data,1) size(data,2) frequencies*phases*wavelengths trials]);
data = mean(data,4);

% create mask of brain
figure
imagesc(squeeze(data(:,:,1))) 
colormap gray
colorbar
axis image
button = 1;
count = 1;
while (button == 1)
    [x(count),y(count),button] = ginput(1);
    hold on
    plot(x(count),y(count),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k')
    count = count+1;
    if button == 3
        break
    end
end
x = round(x);
y = round(y);
bw = poly2mask(x,y,size(data,1),size(data,2));
data = data .* bw;
data_ref = data_ref .* bw;

% crop images
figure
colormap gray
imagesc(squeeze(data(:,:,1))) 
axis image
title('Select crop region')
h = imrect(gca, [20 20 150 150]);
addNewPositionCallback(h,@(p) title(mat2str(p,3)));
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(h,fcn);
position = wait(h);
x(1) = int32(position(1));
x(2) = x(1)+int32(position(3));
y(1) = int32(position(2));
y(2) = y(1)+int32(position(4));
data = data(y(1):y(2),x(1):x(2),:);
data_ref = data_ref(y(1):y(2),x(1):x(2),:);

for w = 1:wavelengths
    count = 1;
    for i = w:phases*wavelengths:frequencies*phases*wavelengths
        M_AC_ref(:,:,count,w) = sqrt(2)./3.*(sqrt((data_ref(:,:,i) - data_ref(:,:,i+wavelengths)).^2 + ...
                                (data_ref(:,:,i+wavelengths) - data_ref(:,:,i+wavelengths*2)).^2 + ...
                                (data_ref(:,:,i+wavelengths*2) - data_ref(:,:,i)).^2));

        M_DC_ref(:,:,count,w) = 1./3.*(data_ref(:,:,i) + data_ref(:,:,i+wavelengths) + data_ref(:,:,i+wavelengths*2));


        M_AC(:,:,count,w) = sqrt(2)./3.*(sqrt((data(:,:,i) - data(:,:,i+wavelengths)).^2 + ...
                                (data(:,:,i+wavelengths) - data(:,:,i+wavelengths*2)).^2 + ...
                                (data(:,:,i+wavelengths*2) - data(:,:,i)).^2));

        M_DC(:,:,count,w) = 1./3.*(data(:,:,i) + data(:,:,i+wavelengths) + data(:,:,i+wavelengths*2));

        count = count+1;
    end

    for i = 1:frequencies
        Rd_dc(:,:,i,w) = M_DC(:,:,i,w)./M_DC_ref(:,:,i,w);
        Rd_ac(:,:,i,w) = M_AC(:,:,i,w)./M_AC_ref(:,:,i,w);
    end
end

%% Use LUT to get reduced scattering and absorption coefficients

kx = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 1 1.2 1.5];
idx1 = 1;   
idx2 = 8;
muaall = 0:0.1:1.3;
musall = 1:3:40;
muasmall = muaall(1):0.01:muaall(end);
mussmall = musall(1):0.3:musall(end);
[reReflect_ac, reReflect_dc] = generateLUT(kx(idx1), muaall, musall);
muainterp = interp1(muaall,reReflect_dc,muasmall);
musinterp = interp1(musall,muainterp',mussmall);
reReflect_dc_interp = musinterp';

[reReflect_ac1, na] = generateLUT(kx(idx2), muaall, musall);
muainterp = interp1(muaall,reReflect_ac1,muasmall);
musinterp = interp1(musall,muainterp',mussmall);
reReflect_ac1_interp = musinterp';

reRef_dc = reshape(reReflect_dc_interp,[size(reReflect_dc_interp,1)*size(reReflect_dc_interp,2),1,1]);
reRef_ac1 = reshape(reReflect_ac1_interp,[size(reReflect_ac1_interp,1)*size(reReflect_ac1_interp,2),1,1]);
mua = repmat(muasmall, [1,length(muasmall)])';
mus = repmat(mussmall, [length(mussmall),1]);
mus = reshape(mus, [size(mus,1)*size(mus,2),1]);
LUTtable = [reRef_dc, reRef_ac1, mua, mus];

for w = 1:wavelengths
    for i = 1:size(Rd_dc,1)
        for j = 1:size(Rd_dc,2)
            if isnan(Rd_dc(i,j,1,w)) == 1 || Rd_dc(i,j,1,w) == Inf
                prop_mua(i,j,w) = 0;
                prop_mus(i,j,w) = 0;
            else
                dist_dc = (LUTtable(:,1) - Rd_dc(i,j,idx1,w)).^2;
                dist_ac = (LUTtable(:,2) - Rd_ac(i,j,idx2,w)).^2;
                dist = (dist_ac + dist_dc);
                prop_mua(i,j,w) = mua(find(dist == min(dist)));
                prop_mus(i,j,w) = mus(find(dist == min(dist)));
                distance(i,j,w) = dist(find(dist == min(dist)));
            end
        end
    end
end

figure
for i = 1:wavelengths
    imagesc(prop_mua(:,:,i))
    axis image
    colorbar
    colormap jet
    caxis([0 1])
    title([num2str(lambda(i)),': \mu_a'])
    axis off
    figure
    imagesc(prop_mus(:,:,i)./10)
    colorbar
    caxis([0.4 2.5])
    axis image
    colormap jet
    title([num2str(lambda(i)),': \mu_s'])
    axis off
end
save([pathname, '\mus.mat'],'prop_mus')
