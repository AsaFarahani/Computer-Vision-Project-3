clc
close all
clear all
%% read images 

I1 = imread("./stereo-pairs/p21.jpg"); %left image
I2 = imread("./stereo-pairs/p22.jpg"); %right image

I1 = rgb2gray(I1);
I2 = rgb2gray(I2);
%% SURF feature detector 
%returns a SURFPoints object, points, containing information about SURF features detected in the 2-D grayscale input image I. 
ptsLeft  = detectSURFFeatures(I1); 
ptsRight = detectSURFFeatures(I2);

%ExtractFeatures: Returns [features, validPoints], extracted feature vectors, also known as
%descriptors, and their corresponding locations.
[featuresLeft,  validPtsLeft]  = extractFeatures(I1,  ptsLeft);
[featuresRight, validPtsRight] = extractFeatures(I2, ptsRight);

%indexPairs the matching feature vectors and return the matching points
%index in featuresLeft, and featuresRight
indexPairs = matchFeatures(featuresLeft, featuresRight);

%Get the valid points using the index.
matchedLeft  = validPtsLeft(indexPairs(:,1));
matchedRight = validPtsRight(indexPairs(:,2));
%% found features in image 1 and image 2 and plot them
features1 = matchedLeft.Location;
features2 = matchedRight.Location;

% plot found pairs 
figure(1);
imshowpair(I1, I2, 'montage');
for i = 1:1:size(features1,1)
    pos1 = features1(i,1:2);
    pos2 = features2(i,1:2);
    pos2(1) = pos2(1)+size(I1, 2);
    x1 = pos1(1); x2 = pos2(1);
    y1 = pos1(2); y2 = pos2(2);
    hold on
    plot([x1, x2], [y1, y2], 'Marker','*');
end
title('data point pairs - SURF features');

%close all % it should be commented 
%% ransac
% my own code of ransac - I have two versions 
F_best = F_ransac(features1,features2)
%F_best = F_ransac1(features1,features2)    
% print rank of matrix F_best
fprintf ('rank of F obtained from ransac is = \n %f \n', rank(F_best)); 

% build in function of matlab - for comparison 
[F_best_matalb,inliersIndex] =estimateFundamentalMatrix(features1,...
    features2,'Method','RANSAC',...
    'NumTrials',10000,'DistanceThreshold',0.01);
F_best_matalb

%% epipoles
e1 = null(F_best); % right null space of F -> epipole 1
e1 = e1./e1(3);
e2 = null(F_best_matalb'); % left null space F -> epipole 2
e2 = e2./e2(3);
%% Homography - sent epipole to infinity
Homography1 = [1 , 0 , 0;...
    -(e1(2)/e1(1)), 1 , 0; ...
    -(1/e1(1)), 0 , 1];
Homography2 = [1 , 0 , 0;...
    -(e2(2)/e2(1)), 1 , 0; ...
    -(1/e2(1)), 0 , 1];
%% apply homography - image 1 
for i = 1:1:size(I1,1)
    for j = 1:1: size(I1,2)
        homography1_applied{i,j} = Homography1 * double([i; j; 1]);
        temp1 = homography1_applied{i,j};
        homography1_applied_norm{i,j} = homography1_applied{i,j}./temp1(3);
    end
end
for i=1:1:size(I1,1)
    for j =1:1:size(I1,2)
        Homo_I1_y(i,j) = homography1_applied_norm{i,j}(1);
        Homo_I1_x(i,j) = homography1_applied_norm{i,j}(2);
    end
end
%% apply homography - image 2
for i = 1:1: size(I2,1)
    for j = 1:1: size(I2,2)
        homography2_applied{i,j} = Homography2 * double([i; j; 1]);
        temp2 = homography2_applied{i,j};
        homography2_applied{i,j} = homography2_applied{i,j}./temp2(3);
    end
end
for i=1:1:size(I2,1)
    for j =1:1: size(I2,2)
        Homo_I2_y(i,j) = homography2_applied{i,j}(1);
        Homo_I2_x(i,j) = homography2_applied{i,j}(2);
    end
end
%% homography - image 1 and 2 together 
%determine the minimum and maximum bounds for the composite image based off
%the warped corners
% img 1 
minY_I1 = min(Homo_I1_y(:));
maxY_I1 = max(Homo_I1_y(:));
minX_I1 = min(Homo_I1_x(:));
maxX_I1 = max(Homo_I1_x(:));
% img2 
minY_I2 = min(Homo_I2_y(:));
maxY_I2 = max(Homo_I2_y(:));
minX_I2 = min(Homo_I2_x(:));
maxX_I2 = max(Homo_I2_x(:));
% overall of img1 and img2 
minx = min(minX_I2,minX_I1);
miny = min(minY_I2,minY_I1);
maxx = max(maxX_I2,maxX_I1);
maxy = max(maxY_I2,maxY_I1);

% use those min and max bounds to define the resolution of the composite image
xRange = minx : maxx; %the range for x pixels
yRange = miny : maxy; %the range for y pixels
[x,y] = meshgrid(xRange,yRange);

% homography inverse 
Homography1_inv = inv(Homography1);
Homography2_inv = inv(Homography2);

warpedHomoScaleFactor = Homography2_inv(1,3) * x + Homography2_inv(2,3) * y + Homography2_inv(3,3);
warpX = (Homography2_inv(1,1) * x + Homography2_inv(2,1) * y + Homography2_inv(3,1)) ./ warpedHomoScaleFactor ;
warpY = (Homography2_inv(1,2) * x + Homography2_inv(2,2) * y + Homography2_inv(3,2)) ./ warpedHomoScaleFactor ;

rwarpedHomoScaleFactor = Homography1_inv(1,3) * x + Homography1_inv(2,3) * y + Homography1_inv(3,3);
rawarpX = (Homography1_inv(1,1) * x + Homography1_inv(2,1) * y + Homography1_inv(3,1)) ./ rwarpedHomoScaleFactor ;
rawarpY = (Homography1_inv(1,2) * x + Homography1_inv(2,2) * y + Homography1_inv(3,2)) ./ rwarpedHomoScaleFactor ;

% interpolate pixel value
first_img  = interp2( im2double(I1), rawarpX, rawarpY, 'cubic') ;
second_img = interp2( im2double(I2), warpX, warpY, 'cubic') ;

% handle NANs 
second_img(isnan(second_img)) = 0 ;
first_img(isnan(first_img)) = 0 ;

% add the blendedLeft and Right halves together while dividing by the
% blendWeight for that pixel.
blendWeight = ~ isnan(second_img) + ~isnan(first_img);
composite = (second_img + first_img) ./ blendWeight;
%{
% show and save results 
imwrite(second_img,'Rect_pair2_2.jpg');
imwrite(first_img,'Rect_pair2_1.jpg');
imwrite(composite,'Rect_pair2_1_share.jpg');
%}
figure, imshow(composite);
title('Alignment by homography- both');
