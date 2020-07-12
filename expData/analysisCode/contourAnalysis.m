%% Contour analysis 

%% Calculate percent change

pathname{1} = 'E:\Scat13_SS64\SFDI\Baseline\mus.mat';
pathname{2} = 'E:\Scat13_SS64\SFDI\24hr\mus_trans.mat';
for i = 1:length(pathname)
    mus(i).prop = load(pathname{i});
    if i == 1
        mus(i).prop.filt = imgaussfilt(mus(i).prop.prop_mus,2);
        rmusSFDI(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    else
        mus(i).prop.prop_mus = fliplr(imrotate(squeeze(mus(i).prop.angio(1,:,:)),270));
        mus(i).prop.filt = imgaussfilt(mus(i).prop.prop_mus,2);
        rmusSFDI(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    end
end

pathname{1} = 'E:\Scat13_SS64\OCT\Baseline\OCTattenuation_trans_SFDI.mat';
pathname{2} = 'E:\Scat13_SS64\OCT\24hr\OCTattenuation_trans.mat';
for i = 1:length(pathname)
    mus(i).prop = load(pathname{i});
    if i == 1
        mus(i).prop.filt = imgaussfilt(fliplr(imrotate((squeeze(mus(i).prop.angio(1,:,:))),270)),2);
        rmusOCT(i).prop = mus(i).prop.filt./mus(1).prop.filt;
    else
        mus(i).prop.prop_mus = fliplr(imrotate(squeeze(mus(i).prop.angio(1,:,:)),270));
        mus(i).prop.filt = imgaussfilt(mus(i).prop.prop_mus,2);
        rmusOCT(i).prop = mus(i).prop.filt./mus(1).prop.filt;
        rmusOCT(i).prop = (rmusOCT(i).prop);
    end
end


figure
subplot(1,2,1)
imagesc(rmusSFDI(2).prop)
axis image
axis off
colormap jet
caxis([0.5 1.5])
subplot(1,2,2)
imagesc(rmusOCT(2).prop)
axis image
axis off
colormap jet
caxis([0 2])
title('Select crop region')
h = imrect(gca, [20 20 150 150]);
addNewPositionCallback(h,@(p) title(mat2str(p,3)));
fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(h,fcn);
position = wait(h);
x(1) = round(position(1));
x(2) = x(1)+round(position(3));
y(1) = round(position(2));
y(2) = y(1)+round(position(4));
rmusOCT(2).prop = rmusOCT(2).prop(y(1):y(2),x(1):x(2),:);
rmusSFDI(2).prop = rmusSFDI(2).prop(y(1):y(2),x(1):x(2),:);


%% SFDI contour

I_SFDI = rmusSFDI(2).prop;

filter = fspecial('average',15);
I_SFDI = filter2(filter,I_SFDI,'same');
figure
imagesc(I_SFDI)
axis image
colormap jet
caxis([0.5 1.5])
prompt = sprintf('Do you want to remove any edge outliers?');
button = questdlg(prompt,'Remove outliers','Yes','No','No');
if strcmpi(button,'Yes')
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
    bw = poly2mask(x,y,size(I_SFDI,1),size(I_SFDI,2));
    I_SFDI = I_SFDI .* bw;
end


I_SFDI(isnan(I_SFDI)) = 0;
I_SFDI(isinf(I_SFDI)) = 0;
figure
imagesc(I_SFDI)
axis image
colormap jet
caxis([0.5 1.5])
title('Draw rectangle')
count = 1;
select = 1;

while select == 1
prompt1 = sprintf('Do you want to select another region with stroke core?');
button1 = questdlg(prompt1,'Select region','Yes','No','No');
    if strcmpi(button1,'Yes')
        change = 1;
        while change == 1
            hold on
            r = drawrectangle;
            mask = createMask(r);
            outline_SFDI{count} = activecontour(I_SFDI, mask, 200, 'Chan-Vese','SmoothFactor',0.2,'ContractionBias',0.5);
            v = visboundaries(outline_SFDI{count},'Color','k');
            prompt2 = sprintf('Do you want to change your selection');
            button2 = questdlg(prompt2,'Change selection','Yes','No','No');
            if strcmpi(button2,'Yes')
                set(r,'Visible','Off')
                outline_SFDI{count} = zeros(size(outline_SFDI{count}));
                set(v,'Visible','Off')
            elseif strcmpi(button2,'No')
                change = 0;          
                count = count+1;
            end
        end
        
    elseif strcmpi(button1,'No')
        select = 2;
    end
end

full_SFDI = zeros(size(outline_SFDI{1}));
for i = 1:length(outline_SFDI)
    full_SFDI = full_SFDI+outline_SFDI{i};
end
idx = find(full_SFDI>0);
full_SFDI(idx) = 1;

figure
imagesc(I_SFDI)
axis image
colormap jet
caxis([0.5 1.5])
hold on
visboundaries(full_SFDI,'Color','k');

%% OCT contour

I_OCT = rmusOCT(2).prop;

filter = fspecial('average',15);
I_OCT = filter2(filter,I_OCT,'same');
figure
imagesc(I_OCT)
axis image
colormap jet
caxis([0 2])
prompt = sprintf('Do you want to remove any edge outliers?');
button = questdlg(prompt,'Remove outliers','Yes','No','No');
if strcmpi(button,'Yes')
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
    bw = poly2mask(x,y,size(I_OCT,1),size(I_OCT,2));
    I_OCT = I_OCT .* bw;
end

I_OCT(isnan(I_OCT)) = 0;
I_OCT(isinf(I_OCT)) = 0;
figure
imagesc(I_OCT)
axis image
colormap jet
caxis([0 2])
title('Draw rectangle')

count = 1;
select = 1;

while select == 1
prompt1 = sprintf('Do you want to select another region with stroke core?');
button1 = questdlg(prompt1,'Select region','Yes','No','No');
    if strcmpi(button1,'Yes')
        change = 1;
        while change == 1
            hold on
            r = drawrectangle;
            mask = createMask(r);
            outline_OCT{count} = activecontour(I_OCT, mask, 200, 'Chan-Vese','SmoothFactor',0.2,'ContractionBias',0.5);
            v = visboundaries(outline_OCT{count},'Color','k');
            prompt2 = sprintf('Do you want to change your selection');
            button2 = questdlg(prompt2,'Change selection','Yes','No','No');
            if strcmpi(button2,'Yes')
                set(r,'Visible','Off')
                outline_OCT{count} = zeros(size(outline_OCT{count}));
                set(v,'Visible','Off')
            elseif strcmpi(button2,'No')
                change = 0;          
                count = count+1;
            end
        end
        
    elseif strcmpi(button1,'No')
        select = 2;
    end
end

full_OCT = zeros(size(outline_OCT{1}));
for i = 1:length(outline_OCT)
    full_OCT = full_OCT+outline_OCT{i};
end
idx = find(full_OCT>0);
full_OCT(idx) = 1;

figure
imagesc(I_OCT)
axis image
colormap jet
caxis([0 2])
hold on
visboundaries(full_OCT,'Color','k');


%% Overlap between outline of SFDI and OCT with Dice's coefficient 
BW1 = full_SFDI;
BW2 = full_OCT;

similarity = dice(BW2,BW1);

%% Save contour and overlap
save(['E:\Scat13_SS64\SFDI\24hr\','overlap.mat'], 'I_SFDI','I_OCT','full_SFDI','full_OCT','similarity')
