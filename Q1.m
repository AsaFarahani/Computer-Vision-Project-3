clc
close all
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load images
I1 = imread("./stereo-pairs/p11.jpg"); %left image
I2 = imread("./stereo-pairs/p12.jpg"); %right image

%I1 = imread("./stereo-pairs/p21.jpg"); %left image
%I2 = imread("./stereo-pairs/p22.jpg"); %right image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data of the 8 pairs

xPos_img1 = [115,295,268,305,270,497,525,682];
yPos_img1 = [152,87,157,194,373,102,185,312];
xPos_img2 = [10,449,172,256,178,540,475,719];
yPos_img2 = [159,78,158,192,364,86,176,311];

%{
xPos_img1 = [211,360,297,575,535,139,657,412];
yPos_img1 = [325,81,229,380,318,416,106,154];
xPos_img2 = [169,347,283,463,488,62,651,400];
yPos_img2 = [325,84,232,381,318,411,102,155];
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot images and data of the 8 pairs
figure(1); subplot(1,2,1); a = imshow(I1); title('image 1');
datacursormode on
for i =1:1:8
    datatip(a,xPos_img1(i),yPos_img1(i));
end
figure(1);subplot(1,2,2);b = imshow(I2); title('image 2');
for i =1:1:8
    datatip(b,xPos_img2(i),yPos_img2(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
features1 = transpose(vertcat(xPos_img1,yPos_img1))
features2 = transpose(vertcat(xPos_img2,yPos_img2));

% plot 8 selected pairs
figure(2);
imshowpair(I1, I2, 'montage');
for i = 1:1:8
    pos1 = features1(i,1:2);
    pos2 = features2(i,1:2);
    pos2(1) = pos2(1)+size(I1, 2);
    x1 = pos1(1); x2 = pos2(1);
    y1 = pos1(2); y2 = pos2(2);
    hold on
    plot([x1, x2], [y1, y2], 'y-','Marker','*')
end
title('data point pairs in image 1 and image 2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data normalization 
mean_x1 = mean(xPos_img1); 
mean_y1 = mean(yPos_img1);
mean_x2 = mean(xPos_img2);
mean_y2 = mean(yPos_img2);

std_x1 = std(features1(:,1));
std_y1 = std(features1(:,2));

std_x2 = std(features2(:,1));
std_y2 = std(features2(:,2));

M1 = [1/(std_x1) 0 0; 0 1/(std_y1) 0; 0 0 1] * [1 0 -mean_x1; 0 1 -mean_y1; 0 0 1];
M2 = [1/(std_x2) 0 0; 0 1/(std_y2) 0; 0 0 1] * [1 0 -mean_x2; 0 1 -mean_y2; 0 0 1];

% apply normalization
for i = 1:1:8
    a = M1 * [xPos_img1(i); yPos_img1(i); 1];
    x1_norm(i) = a(1); y1_norm(i) = a(2);
    b = M2 * [xPos_img2(i); yPos_img2(i); 1];
    x2_norm(i) = b(1); y2_norm(i) = b(2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fundamantal matrix -> called A at this step (Ax = 0)
A = zeros(9,9);
for i = 1:1:8
    A(i,:) = [(x1_norm(i))*(x2_norm(i)) (y1_norm(i))*(x2_norm(i)) (x2_norm(i)) (x1_norm(i))*(y2_norm(i)) ...
               (y1_norm(i))*(y2_norm(i)) (y2_norm(i)) (x1_norm(i)) (y1_norm(i)) 1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply SVD (on A)
[U,S,V] = svd(A);

% find F - the eigenvector that corresponds to the smallest eigenvalue
F_norm_rank3 = [V(1,9),V(2,9),V(3,9);...
    V(4,9), V(5,9),V(6,9);...
    V(7,9), V(8,9),V(9,9)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rank(F_rank3) = 3 , so we need to make it have rank 2
[UF,SF,VF] = svd(F_norm_rank3);
SF(3,3) = 0 ;
F_norm_rank2 = UF*SF*transpose(VF); % This F has rank 2 

% denormalize F
F = transpose(M2) * F_norm_rank2 * M1;
F_final = F/norm(F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find F using Matlab functions
F_test = estimateFundamentalMatrix(features1,features2,'Method','Norm8Point')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epipoles

% version one - first way of finding epipoles
e1 = null(F_final); % right null space of F -> epipole 1
e1 = e1./e1(3);
e2 = null(F_final'); % left null space F -> epipole 2
e2 = e2./e2(3);

% version two - another way of finding epipoles

[u,d] = eigs(transpose(F_final) * F_final);
uu = u(:,3); %Eigenvector associated with smallest eigenvalue
e11 = uu / uu(3); %Epipole projected to image coordinates 
%this is where the other picture is being taken
[u,d] = eigs((F_final) * transpose(F_final));
uu = u(:,3); %Eigenvector associated with smallest eigenvalue
e22 = uu / uu(3); %Epipole projected to image coordinates 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% epipolar lines
figure(4);subplot(1,2,1); imagesc(I1); title('image1');
figure(4);subplot(1,2,2); imagesc(I2); title('image2');
xx = 0:size(I2,2);
for i =1:size(xPos_img1,2)
    % epipolar line 1
    l1 = [xPos_img2(i) yPos_img2(i) 1]*F_final; l1 = l1./l1(3);
    % 2 * x + y - 300 = 0. #line = [2, 1, -300]
    yy = -l1(1)/l1(2)* xx - l1(3)/l1(2);
    figure(4); subplot(1,2,1);
    hold on
    plot(xx,yy,'LineWidth',1);
    % epipolar line 2
    l2 = F_final*[xPos_img1(i); yPos_img1(i); 1]; l2 = l2./l2(3);
    yy = -l2(1)/l2(2)* xx - l2(3)/l2(2);
    figure(4); subplot(1,2,2);
    hold on
    plot(xx,yy,'LineWidth',1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%