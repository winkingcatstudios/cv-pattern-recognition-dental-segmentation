% HW4

% PART i (optional)
%-------------------------------------------------------------------------
% Adaptive histogram equalization

% Original image
[original,map]=imread('teeth_sample.png');
figure,imshow(original)
title('Original image');
w = waitforbuttonpress;
close

converted = double(original);
converted = converted - min(converted(:));
converted = converted./max(converted(:));
% Play with resizing as desired
%converted = imresize(converted,.5);


% Convert and initial histogram
figure
subplot(1,2,1)
imshow(converted)
title('Original image');
subplot(1,2,2)
imhist(converted,64)
title('Original histogram');
w = waitforbuttonpress;
close

% Adaptive hist equalization and new histogram
preprocessed = adapthisteq(converted);
figure
subplot(1,2,1)
imshow(preprocessed)
title('Preprocessed image');
subplot(1,2,2)
imhist(preprocessed,64)
title('Preprocessed histogram');
w = waitforbuttonpress;
close

% New image to use going forward
figure,imshow(preprocessed)
title('Preprocessed image');
w = waitforbuttonpress;
close

% PART ii 
%-------------------------------------------------------------------------
% Extract bounding boxes/tetrahedra for all teeth in the field of view

% PART iia
%-------------------------------------------------------------------------
% Gap valley detection

% Setup for integrated intensity scan, 20 vert boxes
height = size(preprocessed,1);
width = size(preprocessed,2);
n = width / 20; n = floor(n);
boxStart = 1;
boxEnd = n;

% Calculate all vertical box integrated intensities
sumTemp = 0; i=1; j=1;
IntensityMatrixVert = zeros(height,20);
for k = 1:20    
    for i = 1:height
        for j = boxStart:boxEnd
            sumTemp = sumTemp + preprocessed(i,j);
        end
        IntensityMatrixVert(i,k) = sumTemp;
        sumTemp = 0;
    end
    boxStart = boxStart + n; boxEnd = boxEnd + n;
end

% Display integrated intensity vs image, showing first box
figure
subplot(1,2,1);
imshow(preprocessed,'InitialMagnification', 800);
axis on
hold on;
rectangle('Position', [1 1 n height], 'EdgeColor', 'r' );
title('First Vertical Scan (1/20)');
subplot(1,2,2);
plot(IntensityMatrixVert(:,1)), view(90,90);
title('Integrated Intensity (1/20)');
w = waitforbuttonpress;
close

% Find local minima of integrated intensities
LocalMin = islocalmin(IntensityMatrixVert,1,'MinProminence',2);

% Local min mapped to image sized matrix, at half window width between
% first and last window
boxMid = n/2; boxMid = round(boxMid);
LocalMinMapped = zeros(height,width);
for k = 1:20  
    for i = 1:height
        if LocalMin(i,k) == 1
            if k == 1
                LocalMinMapped(i,1) = 1;
            elseif k == 20
                LocalMinMapped(i,width) = 1;
            else
                LocalMinMapped(i,boxMid) = 1;
            end
        end
    end
    boxMid = boxMid + n;
end

% Plot local minima and y-hat 
yHat = height / 2; yHat = round(yHat); % Start with half down the image
yHat = yHat - 40; % Adjust y-hat according to how radiographs are taken
yHatx = [1 width];
yHaty = [yHat yHat];
figure,imshow(preprocessed);
hold on;
line(yHatx,yHaty,'Color','g','LineStyle','--')
hold on;
spy(sparse(LocalMinMapped), 'r*');
title('Local Minima (red) and Y-hat (green)');
w = waitforbuttonpress;
close


% Equations 1-3 in Jain,Chen paper
% eq 3: Pvi(yi) = (1/sqrt(2*pi*sigma))*(e^((-(yi-yhat)^2)/(sigma^2)))
% eq 2: Pvi(Di) = c(1-(Di/(maxDk))
% eq 1: Pvi (Di,yi) = Pvi(yi)*Pvi(Di)
GapValleyProb3 = zeros(height,width); % for eq3
GapValleyProb2 = zeros(height,20); % for eq2
GapValleyProb1 = zeros(height,width); % for eq1
sigma = 1; c = 1;  % Adjust if corrections are needed

% eq 2 (aka: is it dark enough?)
maxDks = []; 
maxDks = max(IntensityMatrixVert,[],1); 
for k = 1:20  
    for i = 1:height
        GapValleyProb2(i,k) = c*(1-(IntensityMatrixVert(i,k)/maxDks(k)));     
    end
end
% Map equation 2 to image size to compare to GapValleyProb3
boxMid = n/2; boxMid = round(boxMid);
GapValleyProb2Mapped = zeros(height,width);
for k = 1:20  
    for i = 1:height
        if k == 1
            GapValleyProb2Mapped(i,1) = GapValleyProb2(i,k);
        elseif k == 20
            GapValleyProb2Mapped(i,width) = GapValleyProb2(i,k);
        else
            GapValleyProb2Mapped(i,boxMid) = GapValleyProb2(i,k);
        end
    end
    boxMid = boxMid + n;
end

% eq 3 (aka: is it close enough to the assumed midpoint of the image?)
for k = 1:width  
    for i = 1:height
        if LocalMinMapped(i,k) == 1
            GapValleyProb3(i,k) = (1/sqrt(2*pi*sigma)*exp((-(i-yHat)^2)/(sigma^2)));
        end
    end
end

% eq 1 (product of above probabilities)
GapX = []; GapY = [];
for k = 1:width  
    for i = 1:height
        GapValleyProb1(i,k) = GapValleyProb3(i,k) * GapValleyProb2Mapped(i,k);
        if GapValleyProb1(i,k) > 0
            GapX = [GapX; k];
            GapY = [GapY; i];
        end
    end
end


% Plot gap valley based on 3 equation probability
fit = 5; % Higher vals = higher degree polynomial = more tight fit
GapValleyPoly = polyfit(GapX, GapY, fit);
GapValleyVal = polyval(GapValleyPoly, GapX);
figure,imshow(preprocessed);
hold on,plot(GapX,GapValleyVal,'r')
title('Gap Valley Detection');
w = waitforbuttonpress;
close

% PART iib
%-------------------------------------------------------------------------
% Tooth isolation

% Same process as iia, but scan horizontally
% Do this twice, 10 for maxillary teeth and 10 for mandibular teeth

% Maxillary 
% Setup for integrated intensity scan, 10 horizontal boxes above Y-hat
heightMaxi = yHat;
width = size(preprocessed,2);
nMaxi = heightMaxi / 10; nMaxi = floor(nMaxi);
boxStart = 1;
boxEnd = nMaxi;

% Calculate all vertical box integrated intensities
sumTemp = 0; i=1; j=1; k=1;
IntensityMatrixMaxi = zeros(10,width);
for k = 1:10   
    for i = 1:width 
        for j = boxStart:boxEnd
            sumTemp = sumTemp + preprocessed(j,i);
        end
        IntensityMatrixMaxi(k,i) = sumTemp;
        sumTemp = 0;
    end
    boxStart = boxStart + nMaxi; boxEnd = boxEnd + nMaxi;
end

% Display integrated intensity vs image, showing first box
lastBoxHeightMaxi = yHat-nMaxi;
figure
subplot(2,1,1);
imshow(preprocessed,'InitialMagnification', 800);
axis on
hold on;
rectangle('Position', [1 lastBoxHeightMaxi width nMaxi], 'EdgeColor', 'r' );
title('Last Maxillary Scan (10/10)');
subplot(2,1,2);
plot(IntensityMatrixMaxi(10,:));
title('Integrated Intensity (10/10)');
w = waitforbuttonpress;
close

% Mandibular
% Setup for integrated intensity scan, 10 horizontal boxes below Y-hat
heightMand = size(preprocessed,1) - yHat;
width = size(preprocessed,2);
nMand = heightMand / 10; nMand = floor(nMand);
boxStart = yHat;
boxEnd = yHat + nMand;

% Calculate all vertical box integrated intensities
sumTemp = 0; i=1; j=1; k=1;
IntensityMatrixMand = zeros(10,width);
for k = 1:10   
    for i = 1:width 
        for j = boxStart:boxEnd
            sumTemp = sumTemp + preprocessed(j,i);
        end
        IntensityMatrixMand(k,i) = sumTemp;
        sumTemp = 0;
    end
    boxStart = boxStart + nMand; boxEnd = boxEnd + nMand;
end

% Display integrated intensity vs image, showing first box
lastBoxHeightMand = yHat+nMand;
figure
subplot(2,1,1);
imshow(preprocessed,'InitialMagnification', 800);
axis on
hold on;
rectangle('Position', [1 lastBoxHeightMand width nMand], 'EdgeColor', 'r' );
title('Second Mandibular Scan (2/10)');
subplot(2,1,2);
plot(IntensityMatrixMand(2,:));
title('Integrated Intensity (2/10)');
w = waitforbuttonpress;
close

% Find local minima of integrated intensities for maxi and mand
LocalMinMaxi = islocalmin(IntensityMatrixMaxi,2,'MinProminence',2);
LocalMinMand = islocalmin(IntensityMatrixMand,2,'MinProminence',2);

% Local min mapped to above image aboce yhat sized matrix, at half window width between
% first and last window
boxMidMaxi = nMaxi/2; boxMidMaxi = round(boxMidMaxi);
LocalMinMappedMaxi = zeros(height,width);
for k = 1:10  
    for i = 1:width
        if LocalMinMaxi(k,i) == 1
            if k == 1
                LocalMinMappedMaxi(1,i) = 1;
            elseif k == 10
                LocalMinMappedMaxi(heightMaxi,i) = 1;
            else
                LocalMinMappedMaxi(boxMidMaxi,i) = 1;
            end
        end
    end
    boxMidMaxi = boxMidMaxi + nMaxi;
end

boxMidMand = yHat + nMand/2; boxMidMand = round(boxMidMand);
LocalMinMappedMand = zeros(height,width);
for k = 1:10  
    for i = 1:width
        if LocalMinMand(k,i) == 1
            if k == 1
                LocalMinMappedMand(yHat,i) = 1;
            elseif k == 10
                LocalMinMappedMand(height,i) = 1;
            else
                LocalMinMappedMand(boxMidMand,i) = 1;
            end
        end
    end
    boxMidMand = boxMidMand + nMand;
end

% Plot local horizontal minima
figure,imshow(preprocessed);
hold on;
spy(sparse(LocalMinMappedMaxi), 'r*');
hold on;
spy(sparse(LocalMinMappedMand), 'r*');
title('Local Minima (red) and X-hats (green)');
w = waitforbuttonpress;
close

% Manually set X-hats, as allowed by assignment
% This could be done programmatically by taking a mid scan (5/10)
% integrated intesities of Maxi and Mandi and finding the dips to give 
% a rough idea of the X-hat since this is just an estimate
% numXhats = 1 per gap betweeen teeth. Adjust based on number of expected teeth, depends on the type of radiograph
numXhats = 9;  % 11 expected teeth in image, 9 gaps expected
figure,imshow(preprocessed)
title('Select X-hats (gaps between visible teeth)');
[x,y]=ginput(numXhats);
close

% Plot local Maxillary minima and X-hats
figure,imshow(preprocessed);
hold on;
xline(x(1),'Color','g','LineStyle','--')
hold on;
xline(x(2),'Color','g','LineStyle','--')
hold on;
xline(x(3),'Color','g','LineStyle','--')
hold on;
xline(x(4),'Color','g','LineStyle','--')
hold on;
spy(sparse(LocalMinMappedMaxi), 'r*');
title('Local Maxillary Minima (red) and X-hats (green)');
w = waitforbuttonpress;
close

% Plot local Mandibular minima and X-hats
figure,imshow(preprocessed);
hold on;
xline(x(5),'Color','g','LineStyle','--')
hold on;
xline(x(6),'Color','g','LineStyle','--')
hold on;
xline(x(7),'Color','g','LineStyle','--')
hold on;
xline(x(8),'Color','g','LineStyle','--')
hold on;
xline(x(9),'Color','g','LineStyle','--')
hold on;
spy(sparse(LocalMinMappedMand), 'r*');
title('Local Mandibular Minima (red) and X-hats (green)');
w = waitforbuttonpress;
close

% Equations for Maxi
% Equations 1-3 in Jain,Chen paper
% eq 3: Pvi(yi) = (1/sqrt(2*pi*sigma))*(e^((-(yi-yhat)^2)/(sigma^2)))
% eq 2: Pvi(Di) = c(1-(Di/(maxDk))
% eq 1: Pvi (Di,yi) = Pvi(yi)*Pvi(Di)
GapValleyProb3Maxi = zeros(height,width); % for eq3
GapValleyProb2Maxi = zeros(height,20); % for eq2
GapValleyProb1Maxi = zeros(height,width); % for eq1
sigma = 1; c = 1;  % Adjust if corrections are needed

% eq 2 (aka: is it dark enough?)
maxDksMaxi = []; 
maxDksMaxi = max(IntensityMatrixMaxi,[],2); 
for k = 1:10  
    for i = 1:width
        GapValleyProb2Maxi(k,i) = c*(1-(IntensityMatrixMaxi(k,i)/maxDksMaxi(k)));     
    end
end
% Map equation 2 to image size to compare to GapValleyProb3Maxi
boxMidMaxi = nMaxi/2; boxMidMaxi = round(boxMidMaxi);
GapValleyProb2MaxiMapped = zeros(height,width);
for k = 1:10  
    for i = 1:width
        if k == 1
            GapValleyProb2MaxiMapped(1,i) = GapValleyProb2Maxi(k,i);
        elseif k == 10
            GapValleyProb2MaxiMapped(heightMaxi,i) = GapValleyProb2Maxi(k,i);
        else
            GapValleyProb2MaxiMapped(boxMidMaxi,i) = GapValleyProb2Maxi(k,i);
        end
    end
    boxMidMaxi = boxMidMaxi + nMaxi;
end

% eq 3 (aka: is it close enough to the assumed midpoint of the image?)
% Ajusted for multiple x-hats
dev = round(x(1)/5); % acceptable X-hat deviation, adjust based on radiograph
for k = 1:height  
    for i = 1:width
        if LocalMinMappedMaxi(k,i) == 1            
            if i > x(1)-dev && i < x(1)+dev % x-hat 1
                GapValleyProb3Maxi(k,i) = LocalMinMappedMaxi(k,i);
            elseif i > x(2)-dev && i < x(2)+dev % x-hat 2
                GapValleyProb3Maxi(k,i) = LocalMinMappedMaxi(k,i);
            elseif i > x(3)-dev && i < x(3)+dev % x-hat 3
                GapValleyProb3Maxi(k,i) = LocalMinMappedMaxi(k,i);
            elseif i > x(4)-dev && i < x(4)+dev % x-hat 4
                GapValleyProb3Maxi(k,i) = LocalMinMappedMaxi(k,i);
            end    
        end
    end
end

% eq 1 (product of above probabilities)
GapXMaxi1 = x(1); GapYMaxi1 = 1; % x-hat 1
GapXMaxi2 = x(2); GapYMaxi2 = 1; % x-hat 2
GapXMaxi3 = x(3); GapYMaxi3 = 1; % x-hat 3
GapXMaxi4 = x(4); GapYMaxi4 = 1; % x-hat 4
for k = 1:height  
    for i = 1:width
        GapValleyProb1Maxi(k,i) = GapValleyProb3Maxi(k,i) * GapValleyProb2MaxiMapped(k,i);
        if GapValleyProb1Maxi(k,i) > 0
            if i > x(1)-dev && i < x(1)+dev % x-hat 1
                GapXMaxi1 = [GapXMaxi1; i];
                GapYMaxi1 = [GapYMaxi1; k];
            elseif i > x(2)-dev && i < x(2)+dev % x-hat 2
                GapXMaxi2 = [GapXMaxi2; i];
                GapYMaxi2 = [GapYMaxi2; k];
            elseif i > x(3)-dev && i < x(3)+dev % x-hat 3
                GapXMaxi3 = [GapXMaxi3; i];
                GapYMaxi3 = [GapYMaxi3; k];
            elseif i > x(4)-dev && i < x(4)+dev % x-hat 4
                GapXMaxi4 = [GapXMaxi4; i];
                GapYMaxi4 = [GapYMaxi4; k];
            end   
        end
    end
end

% Extend to gap valley from part iia, this is just to make it look nicer so
% bourders mostly intersect. This is manual based on gap valley Y vals
GapXMaxi1 = [GapXMaxi1; GapXMaxi1(length(GapXMaxi1))];
GapYMaxi1 = [GapYMaxi1; GapY(4)];
GapXMaxi2 = [GapXMaxi2; GapXMaxi2(length(GapXMaxi2))];
GapYMaxi2 = [GapYMaxi2; GapY(7)];
GapXMaxi3 = [GapXMaxi3; GapXMaxi3(length(GapXMaxi3))];
GapYMaxi3 = [GapYMaxi3; GapY(10)];
GapXMaxi4 = [GapXMaxi4; GapXMaxi4(length(GapXMaxi4))];
GapYMaxi4 = [GapYMaxi4; GapY(16)];

% Plot gap valley based on 3 equation probability
figure,imshow(preprocessed);
hold on,line(GapXMaxi1,GapYMaxi1,'Color','r')
hold on,line(GapXMaxi2,GapYMaxi2,'Color','r')
hold on,line(GapXMaxi3,GapYMaxi3,'Color','r')
hold on,line(GapXMaxi4,GapYMaxi4,'Color','r')
hold on,plot(GapX,GapValleyVal,'r')
title('Gap Valley Detection');
w = waitforbuttonpress;
close

% Equations for Mand
% Equations 1-3 in Jain,Chen paper
% eq 3: Pvi(yi) = (1/sqrt(2*pi*sigma))*(e^((-(yi-yhat)^2)/(sigma^2)))
% eq 2: Pvi(Di) = c(1-(Di/(maxDk))
% eq 1: Pvi (Di,yi) = Pvi(yi)*Pvi(Di)
GapValleyProb3Mand = zeros(height,width); % for eq3
GapValleyProb2Mand = zeros(height,20); % for eq2
GapValleyProb1Mand = zeros(height,width); % for eq1
sigma = 1; c = 1;  % Adjust if corrections are needed

% eq 2 (aka: is it dark enough?)
maxDksMand = []; 
maxDksMand = max(IntensityMatrixMand,[],2); 
for k = 1:10  
    for i = 1:width
        GapValleyProb2Mand(k,i) = c*(1-(IntensityMatrixMand(k,i)/maxDksMand(k)));     
    end
end
% Map equation 2 to image size to compare to GapValleyProb3Mand
boxMidMand = yHat + nMand/2; boxMidMand = round(boxMidMand);
GapValleyProb2MandMapped = zeros(height,width);
for k = 1:10  
    for i = 1:width
        if k == 1
            GapValleyProb2MandMapped(yHat,i) = GapValleyProb2Mand(k,i);
        elseif k == 10
            GapValleyProb2MandMapped(height,i) = GapValleyProb2Mand(k,i);
        else
            GapValleyProb2MandMapped(boxMidMand,i) = GapValleyProb2Mand(k,i);
        end
    end
    boxMidMand = boxMidMand + nMand;
end

%%%%%%%%%%
% boxMidMand = yHat + nMand/2; boxMidMand = round(boxMidMand);
% LocalMinMappedMand = zeros(height,width);
% for k = 1:10  
%     for i = 1:width
%         if LocalMinMand(k,i) == 1
%             if k == 1
%                 LocalMinMappedMand(yHat,i) = 1;
%             elseif k == 10
%                 LocalMinMappedMand(height,i) = 1;
%             else
%                 LocalMinMappedMand(boxMidMand,i) = 1;
%             end
%         end
%     end
%     boxMidMand = boxMidMand + nMand;
% end
%%%%%%%%%%%%%%

% eq 3 (aka: is it close enough to the assumed midpoint of the image?)
% Ajusted for multiple x-hats
dev = round(x(1)/5); % acceptable X-hat deviation, adjust based on radiograph
for k = 1:height  
    for i = 1:width
        if LocalMinMappedMand(k,i) == 1            
            if i > x(5)-dev && i < x(5)+dev % x-hat 5
                GapValleyProb3Mand(k,i) = LocalMinMappedMand(k,i);
            elseif i > x(6)-dev && i < x(6)+dev % x-hat 6
                GapValleyProb3Mand(k,i) = LocalMinMappedMand(k,i);
            elseif i > x(7)-dev && i < x(7)+dev % x-hat 7
                GapValleyProb3Mand(k,i) = LocalMinMappedMand(k,i);
            elseif i > x(8)-dev && i < x(8)+dev % x-hat 8
                GapValleyProb3Mand(k,i) = LocalMinMappedMand(k,i);
            elseif i > x(9)-dev && i < x(9)+dev % x-hat 9
                GapValleyProb3Mand(k,i) = LocalMinMappedMand(k,i);
            end    
        end
    end
end

% eq 1 (product of above probabilities)
GapXMand5 = x(5); GapYMand5 = yHat; % x-hat 1
GapXMand6 = x(6); GapYMand6 = yHat; % x-hat 2
GapXMand7 = x(7); GapYMand7 = yHat; % x-hat 3
GapXMand8 = x(8); GapYMand8 = yHat; % x-hat 4
GapXMand9 = x(9); GapYMand9 = yHat; % x-hat 5
for k = 1:height  
    for i = 1:width
        GapValleyProb1Mand(k,i) = GapValleyProb3Mand(k,i) * GapValleyProb2MandMapped(k,i);
        if GapValleyProb1Mand(k,i) > 0
            if i > x(5)-dev && i < x(5)+dev % x-hat 1
                GapXMand5 = [GapXMand5; i];
                GapYMand5 = [GapYMand5; k];
            elseif i > x(6)-dev && i < x(6)+dev % x-hat 2
                GapXMand6 = [GapXMand6; i];
                GapYMand6 = [GapYMand6; k];
            elseif i > x(7)-dev && i < x(7)+dev % x-hat 3
                GapXMand7 = [GapXMand7; i];
                GapYMand7 = [GapYMand7; k];
            elseif i > x(8)-dev && i < x(8)+dev % x-hat 4
                GapXMand8 = [GapXMand8; i];
                GapYMand8 = [GapYMand8; k];
            elseif i > x(9)-dev && i < x(9)+dev % x-hat 4
                GapXMand9 = [GapXMand9; i];
                GapYMand9 = [GapYMand9; k];
            end   
        end
    end
end

% % Extend to gap valley from part iia, this is just to make it look nicer so
% % bourders mostly intersect. This is manual based on gap valley Y vals
GapXMand5 = [GapXMand5; GapXMand5(length(GapXMand5))];
GapYMand5 = [GapYMand5; height];
GapXMand6 = [GapXMand6; GapXMand6(length(GapXMand6))];
GapYMand6 = [GapYMand6; height];
GapXMand7 = [GapXMand7; GapXMand7(length(GapXMand7))];
GapYMand7 = [GapYMand7; height];
GapXMand8 = [GapXMand8; GapXMand8(length(GapXMand8))];
GapYMand8 = [GapYMand8; height];
GapXMand9 = [GapXMand9; GapXMand9(length(GapXMand9))];
GapYMand9 = [GapYMand9; height];

% Plot gap valley based on 3 equation probability
figure,imshow(preprocessed);
hold on,line(GapXMand5,GapYMand5,'Color','r')
hold on,line(GapXMand6,GapYMand6,'Color','r')
hold on,line(GapXMand7,GapYMand7,'Color','r')
hold on,line(GapXMand8,GapYMand8,'Color','r')
hold on,line(GapXMand9,GapYMand9,'Color','r')
hold on,plot(GapX,GapValleyVal,'r')
title('Gap Valley Detection');
w = waitforbuttonpress;
close


% Finished Bounding Boxes
figure,imshow(preprocessed);
hold on,line(GapXMaxi1,GapYMaxi1,'Color','r')
hold on,line(GapXMaxi2,GapYMaxi2,'Color','r')
hold on,line(GapXMaxi3,GapYMaxi3,'Color','r')
hold on,line(GapXMaxi4,GapYMaxi4,'Color','r')
hold on,line(GapXMand5,GapYMand5,'Color','r')
hold on,line(GapXMand6,GapYMand6,'Color','r')
hold on,line(GapXMand7,GapYMand7,'Color','r')
hold on,line(GapXMand8,GapYMand8,'Color','r')
hold on,line(GapXMand9,GapYMand9,'Color','r')
hold on,plot(GapX,GapValleyVal,'r')
title('Teeth Bounding Boxes');
w = waitforbuttonpress;
close


