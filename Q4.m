clc
close all
clear all
%% read images 

I1 = imread("./stereo-pairs/p21.jpg"); %left image
I2 = imread("./stereo-pairs/p22.jpg"); %right image
%{
I1 = imread("./stereo-pairs/p11.jpg"); %left image
I2 = imread("./stereo-pairs/p12.jpg"); %right image
%}
I1 = rgb2gray(I1);
I2 = rgb2gray(I2);
%% intrinsic matrix 
f = 32
s = 31.2*((10)^(-3))
px = round(size(I1,1)/2);
py = round(size(I1,2)/2);

K = [f/s 0 px; 0 f/s py; 0  0 1];
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
%% ransac - find F 
% my own code of ransac - I have two versions 
F = F_ransac(features1,features2)
%% Essential matrix
E =  K' * F * K;
%% Decompose E to find rotation and translation
[U,S,V]= svd(E);
e = (S(1,1) + S(2,2)) / 2;
S(1,1) = e;
S(2,2) = e;
S(3,3) = 0;
E = U * S * V';

[U, ~, V] = svd(E);
W = [0 -1 0; 1 0 0; 0 0 1];
Z = [0 1 0; -1 0 0; 0 0 0];

% Two possible rotation matrices
R1 = U * W * V';
R2 = U * W' * V';

% Force rotations to have the needed properties
if det(R1) < 0
    R1 = -R1;
end
if det(R2) < 0
    R2 = -R2;
end
% Two possible translation vectors
S1 =  U(:,3);
S2 = -U(:,3);

%% Which P matrix is the best version? 
P(:,:,1) =  K*[R1, S1];
P(:,:,2) =  K*[R1, S2];
P(:,:,3) =  K*[R2, S1];
P(:,:,4) =  K*[R2, S2];
%{
xPos_img1 = [270,279,115,295,268,305,270,497,525,682];
yPos_img1 = [375,352,152,87,157,194,373,102,185,312];
xPos_img2 = [178,233,10,449,172,256,178,540,475,719];
yPos_img2 = [364,343,159,78,158,192,364,86,176,311];
%}
xPos_img1 = [211,297,575,535,139];
yPos_img1 = [325,229,380,318,416];
xPos_img2 = [169,283,463,488,62];
yPos_img2 = [325,232,381,318,411];

clear features1
clear features2
features1 = transpose(vertcat(xPos_img1,yPos_img1))
features2 = transpose(vertcat(xPos_img2,yPos_img2));

%%
P1 = K*horzcat([1,0,0;0,1,0;0,0,1], [0;0;0])
P2 = P

X1 = [features1 ones(size(features1,1),1)];   % homogeneous coordinate
X2 = [features2 ones(size(features2,1),1)];

x11 = X1';
x22 = X2';
for m = 1:1:4
    for i=1:size(x11,2)
        %Select point
        sx1 = x11(:,i);
        sx2 = x22(:,i);
             A = [sx1(1,1).*P1(3,:) - P1(1,:);...
             sx1(2,1).*P1(3,:) - P1(2,:);...
             sx2(1,1).*P2(3,:,m) - P2(1,:,m);...
             sx2(2,1).*P2(3,:,m) - P2(2,:,m);]

        % SDV(A)
        [U,D,V] = svd(A);   
       
        %Point in 3D space is the last column of V
        X_temp = V(:,4);
        points_in_3D_temp = X_temp ./ repmat(X_temp(4,1),4,1);
   
        Xworld(:,i,m) = points_in_3D_temp;
    end
end

% find the best P based on the position of cameras 
for  k = 1:1:4
    figure(k),scatter3(Xworld(1,:,k),Xworld(2,:,k),Xworld(3,:,k));
    % find camera position
    NP1 = null(P1);
    cameraCenter1 = NP1(1:3)/NP1(4)
    
    NP2 = null(P2(:,:,k));
    cameraCenter2 = NP2(1:3)/NP2(4)
    hold on
    scatter3(cameraCenter1(1),cameraCenter1(2),cameraCenter1(3),'filled');
    text(cameraCenter1(1),cameraCenter1(2),cameraCenter1(3),'cameraCenter1')
    scatter3(cameraCenter2(1),cameraCenter2(2),cameraCenter2(3),'filled');
    text(cameraCenter2(1),cameraCenter2(2),cameraCenter2(3),'cameraCenter2')
    hold off
end

% create images + depth values on them
for k =1 :1:4
    figure(4+k)
    imshow(I1) 
    hold on
    for i=1:size(x11,2)
        txt = num2str(Xworld(3,i,k));
        text(X1(i,1),X1(i,2),sin(pi),txt, 'color', 'r')
    end
end