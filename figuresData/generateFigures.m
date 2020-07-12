%% Stroke core revealed by tissue scattering using spatial frequency domain imaging
% Smrithi Sunil1,*, Sefik Evren Erdener1,2, Xiaojun Cheng1, Sreekanth Kura1, Jianbo Tang1, 
% John Jiang1, Kavon Karrobi1, K?v?lc?m K?l?ç1, Darren Roblyer1,3, and David A. Boas1,3
%
% 1Department of Biomedical Engineering, Boston University, Boston, MA 02215, USA
% 2Institute of Neurological Sciences and Psychiatry, Hacettepe University, Ankara, Turkey
% 3Department of Electrical and Computer Engineering, Boston University, Boston, MA 02215, USA
% 
% *Corresponding Author: Smrithi Sunil, E-mail: ssunil@bu.edu

%% Select figures folder
figpath = pwd; %uigetdir('Select figuresData folder');

%% Figure 1C

path = [figpath,'\','figure1C'];
load([path,'\','OCTangio.mat']);
figure('Name','Figure: 1C_1','NumberTitle','off');
imagesc(yCoor,xCoor,xyAG);  
caxis([-1.5 1.5])
colormap('gray')
xlabel( 'X [um]'); ylabel('Y [um]')
axis equal
axis tight

% load([path,'\','OCTintensity.mat']);
load([path,'\','depth.mat']);
figure('Name','Figure: 1C_2','NumberTitle','off');
imshow(depth,[-6 0]); 
colormap gray  

load([path,'\','intDepth.mat']);
figure('Name','Figure: 1C_3','NumberTitle','off');
% z=90:240; % z-level for profiling, 150px for 5x, 50px for 10x
% z=z';
% Z_mm=linspace(0,3.3*length(z)*1e-3,length(z))';
% pos1 = [390 155 25 25];
% int_ROI1 = int(pos1(2):1:pos1(2)+pos1(4),pos1(1):1:pos1(1)+pos1(3),:);
% meanint_ROI1 = mean(mean(int_ROI1,1),2);
% meanint_ROI1 = squeeze(meanint_ROI1(1,1,z));
% fit1 = polyfit(Z_mm,squeeze((meanint_ROI1)),1);
% f1 = polyval(fit1,Z_mm);
plot(Z_mm,f1,'b-')
hold on
% pos2 = [390 200 25 25];
% int_ROI1 = int(pos2(2):1:pos2(2)+pos2(4),pos2(1):1:pos2(1)+pos2(3),:);
% meanint_ROI2 = mean(mean(int_ROI1,1),2);
% meanint_ROI2 = squeeze(meanint_ROI2(1,1,z));
% fit2 = polyfit(Z_mm,squeeze((meanint_ROI2)),1);
% f2 = polyval(fit2,Z_mm);
plot(Z_mm,f2,'g-')
plot(Z_mm,meanint_ROI2,'g.')
plot(Z_mm,meanint_ROI1,'b.')
legend('Within stroke core','Outside stroke core')

load([path,'\','OCTattenuation.mat']);
figure('Name','Figure: 1C_4','NumberTitle','off');
imshow(K, 'DisplayRange', [-7 0])
colormap gray
colorbar 

clearvars -except figpath

%% Figure 2C

path = [figpath,'\','figure2C'];
pathname{1} = [path,'\','prestroke_mus.mat'];
pathname{2} = [path,'\','1hr_mus.mat'];
pathname{3} = [path,'\','2hr_mus.mat'];
pathname{4} = [path,'\','4hr_mus.mat'];
figure('Name','Figure: 2C_1','NumberTitle','off');
for i = 1:length(pathname)
    mus(i).prop = load(pathname{i});
    subplot(1,length(pathname),i)
    imagesc(mus(i).prop.prop_mus./10)
    axis image
    axis off
    colormap jet
    caxis([0.5 2.5])
end

pathname{1} = [path,'\','prestroke_mus.mat'];
pathname{2} = [path,'\','1hr_mus_trans.mat'];
pathname{3} = [path,'\','2hr_mus_trans.mat'];
pathname{4} = [path,'\','4hr_mus_trans.mat'];
figure('Name','Figure: 2C_2','NumberTitle','off');
for i = 1:length(pathname)
    mus(i).prop = load(pathname{i});
    if i == 1
        mus(i).prop.filt = imgaussfilt(mus(i).prop.prop_mus,2);
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    else
        mus(i).prop.prop_mus = fliplr(imrotate(squeeze(mus(i).prop.angio(1,:,:)),270));
        mus(i).prop.filt = imgaussfilt(mus(i).prop.prop_mus,2);
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    end
    subplot(1,length(pathname),i)

    imagesc(rmus(i).prop)
    axis image
    axis off
    colormap jet
    caxis([0.5 1.5])
end
clearvars -except figpath

%% Figure3

% Figure3A
path = [figpath,'\','figure3A'];
load([path,'\','72hr_OCTattenuation.mat']);
figure('Name','Figure: 3A_1','NumberTitle','off');
imshow(K, 'DisplayRange', [-7 0])
colormap gray
colorbar 

pathname{1} = [path,'\','pre_OCTattenuation_trans_SFDI.mat'];
pathname{2} = [path,'\','72hr_OCTattenuation_trans.mat'];
figure('Name','Figure: 3A_2','NumberTitle','off');
for i = 1:length(pathname)
    mus(i).prop = load(pathname{i});
    if i == 1   
        mus(i).prop.filt = fliplr(imrotate((imgaussfilt(squeeze(mus(i).prop.angio(1,:,:)),2)),270));
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    else
        mus(i).prop.prop_mus = squeeze(mus(i).prop.angio(1,:,:));
        mus(i).prop.filt = fliplr(imrotate((imgaussfilt(mus(i).prop.prop_mus,2)),270));
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    end
    imagesc(rmus(i).prop)
    axis image
    axis off
    colormap jet
    caxis([0 2])
end
clearvars -except figpath

% Figure3C
figure('Name','Figure: 3C','NumberTitle','off');
overlap = [0.9136 0.533 0.7469 0.6587 0.7809];
bar(overlap)
set(gca,'Fontsize',14)
xticklabels = {'\begin{tabular}{c} 1 \\ (24hr) \end{tabular}',...
    '\begin{tabular}{c} 2 \\ (24hr) \end{tabular}',...
    '\begin{tabular}{c} 3 \\ (72hr) \end{tabular}',...
    '\begin{tabular}{c} 4 \\ (72hr) \end{tabular}',...
    '\begin{tabular}{c} 5 \\ (72hr) \end{tabular}'};
set(gca,'xtick', 1:5, 'XTickLabel', xticklabels, 'TickLabelInterpreter', 'latex');
clearvars -except figpath

%% Figure4A

path = [figpath,'\','figure4A'];
pathname{1} = [path,'\','prestroke_mus.mat'];
pathname{2} = [path,'\','1hr_mus_trans.mat'];
pathname{3} = [path,'\','2hr_mus_trans.mat'];
pathname{4} = [path,'\','4hr_mus_trans.mat'];
pathname{5} = [path,'\','24hr_mus_trans.mat'];
pathname{6} = [path,'\','72hr_mus_trans.mat'];

figure('Name','Figure: 4A_1','NumberTitle','off');
for i = 1:length(pathname)
    mus(i).prop = load(pathname{i});
    if i == 1
        mus(i).prop.filt = imgaussfilt(mus(i).prop.prop_mus,2);
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    else
        mus(i).prop.prop_mus = fliplr(imrotate(squeeze(mus(i).prop.angio(1,:,:)),270));
        mus(i).prop.filt = imgaussfilt(mus(i).prop.prop_mus,2);
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    end
    subplot(1,length(pathname),i)

    imagesc(rmus(i).prop)
    axis image
    axis off
    colormap jet
    caxis([0.5 1.5])
end
clear pathname 

pathname{1} = [path,'\','prestroke_OCTattenuation_trans.mat'];
pathname{2} = [path,'\','1hr_OCTattenuation_trans.mat'];
pathname{3} = [path,'\','2hr_OCTattenuation_trans.mat'];
pathname{4} = [path,'\','4hr_OCTattenuation_trans.mat'];
pathname{5} = [path,'\','24hr_OCTattenuation_trans.mat'];
pathname{6} = [path,'\','72hr_OCTattenuation_trans.mat'];
figure('Name','Figure: 4A_2','NumberTitle','off');
for i = 1:length(pathname)
    mus(i).prop = load(pathname{i});
    
    if i == 1   
        mus(i).prop.filt = fliplr(imrotate((imgaussfilt(squeeze(mus(i).prop.angio(1,:,:)),2)),270));
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    else
        mus(i).prop.prop_mus = squeeze(mus(i).prop.angio(1,:,:));
        mus(i).prop.filt = fliplr(imrotate((imgaussfilt(mus(i).prop.prop_mus,2)),270));
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    end

    subplot(1,length(pathname),i)
    imagesc(rmus(i).prop)
    axis image
    axis off
    colormap jet
    caxis([0 2])
end
clearvars -except figpath

%% Figure5

path = [figpath,'\','figure5'];
% Figure 5A
pathname{1} = [path,'\','prestroke_mus.mat'];
pathname{2} = [path,'\','2hr_mus_trans.mat'];
figure('Name','Figure: 5A','NumberTitle','off');
for i = 1:length(pathname)
    mus(i).prop = load(pathname{i});
    if i == 1
        mus(i).prop.filt = imgaussfilt(mus(i).prop.prop_mus,2);
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    else
        mus(i).prop.prop_mus = fliplr(imrotate(squeeze(mus(i).prop.angio(1,:,:)),270));
        mus(i).prop.filt = imgaussfilt(mus(i).prop.prop_mus,2);
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    end
end
subplot(1,2,1)
imagesc(rmus(i).prop)
axis image
axis off
colormap jet
caxis([0.5 1.5])
clear pathname 
pathname{1} = [path,'\','prestroke_OCTattenuation_trans.mat'];
pathname{2} = [path,'\','2hr_OCTattenuation_trans.mat'];
for i = 1:length(pathname)
    mus(i).prop = load(pathname{i});
    if i == 1   
        mus(i).prop.filt = fliplr(imrotate((imgaussfilt(squeeze(mus(i).prop.angio(1,:,:)),2)),270));
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    else
        mus(i).prop.prop_mus = squeeze(mus(i).prop.angio(1,:,:));
        mus(i).prop.filt = fliplr(imrotate((imgaussfilt(mus(i).prop.prop_mus,2)),270));
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    end
end
subplot(1,2,2)
imagesc(rmus(i).prop)
axis image
axis off
colormap jet
caxis([0 2])
    
% Figure 5B
load([path,'\','overlap.mat'])
figure('Name','Figure: 5B_1','NumberTitle','off');
subplot(1,2,1)
imagesc(I_SFDI)
hold on
visboundaries(full_SFDI,'Color','k');
caxis([0.5 1.5])
axis image
axis off
colormap jet
subplot(1,2,2)
imagesc(I_OCT)
hold on
visboundaries(full_OCT,'Color','k');
caxis([0 2])
axis image
axis off
colormap jet

figure('Name','Figure: 5B_2','NumberTitle','off');
subplot(1,2,1)
imagesc(full_SFDI)
axis image
axis off
colormap gray
subplot(1,2,2)
imagesc(full_OCT)
axis image
axis off
colormap gray

% Figure 5C
MeanSimilarity = [0.79 0.82 0.68 0.35 0.75];
StdDev = [0.09 0.07 0.13 0.21 0.02];
figure('Name','Figure: 5C','NumberTitle','off');
bar(1:5,MeanSimilarity)
hold on
errorbar(1:5,MeanSimilarity,StdDev,'k','Linestyle','none')
ylim([0 1])
xticklabels = {'\begin{tabular}{c} 1 hour \\ (n=7) \end{tabular}',...
    '\begin{tabular}{c} 2 hours \\ (n=7) \end{tabular}',...
    '\begin{tabular}{c} 4 hours \\ (n=7) \end{tabular}',...
    '\begin{tabular}{c} 24 hours \\ (n=7) \end{tabular}',...
    '\begin{tabular}{c} 72 hours \\ (n=3) \end{tabular}'};
set(gca,'xtick', 1:5, 'XTickLabel', xticklabels, 'TickLabelInterpreter', 'latex','Fontsize',14);
clearvars -except figpath

%% Figure6

path = [figpath,'\','figure6'];
% Figure 6A
pathname{1} = [path,'\','prestroke_mus.mat'];
pathname{2} = [path,'\','24hr_mus_trans.mat'];
load([path,'\','contour.mat']);
figure('Name','Figure: 6A','NumberTitle','off');
for i = 1:length(pathname)
    mus(i).prop = load(pathname{i});
    if i == 1
        mus(i).prop.filt = imgaussfilt(mus(i).prop.prop_mus,2);
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    else
        mus(i).prop.prop_mus = fliplr(imrotate(squeeze(mus(i).prop.angio(1,:,:)),270));
        mus(i).prop.filt = imgaussfilt(mus(i).prop.prop_mus,2);
        rmus(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    end
end
subplot(1,2,1)
imagesc(rmus(i).prop)
axis image
axis off
colormap jet
caxis([0.5 1.5])
subplot(1,2,2)
imagesc(I)
hold on
visboundaries(full,'Color','k');
caxis([0.5 1.5])
axis image
axis off
colormap jet

% Figure 6C
overlap = [0.764 0.6745 0.7068 0.8666 0.6649];
figure('Name','Figure: 6C','NumberTitle','off');
bar(overlap)
set(gca,'Fontsize',14)
xticklabels = {'\begin{tabular}{c} 1 \\ (24hr) \end{tabular}',...
    '\begin{tabular}{c} 2 \\ (24hr) \end{tabular}',...
    '\begin{tabular}{c} 3 \\ (72hr) \end{tabular}',...
    '\begin{tabular}{c} 4 \\ (72hr) \end{tabular}',...
    '\begin{tabular}{c} 5 \\ (72hr) \end{tabular}'};
set(gca,'xtick', 1:5, 'XTickLabel', xticklabels, 'TickLabelInterpreter', 'latex');
clearvars -except figpath
