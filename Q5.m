clc
close all;
clear all;
%% load data
%{
img_1 = imread("./Part 3 - rectified_images/rectified_images_2/Rect_pair2_1.jpg"); %left image
img_2 = imread("./Part 3 - rectified_images/rectified_images_2/Rect_pair2_2.jpg"); %right image
%}
img_1 = imread("./Part 3 - rectified_images/rectified_images_1/Rect_pair_1.jpg"); %left image
img_2 = imread("./Part 3 - rectified_images/rectified_images_1/Rect_pair_2.jpg"); %right image
%% Find the Disparity Map

% define window size
window_size =3;
half_window_size=(window_size-1)/2;

% d - to check
d_min = 0; d_max = 100;

disparity_map = zeros(size(img_1,1),size(img_1,2));

sigma = 6
[map1, map2] = meshgrid(-half_window_size:half_window_size, -half_window_size:half_window_size);
gauss = exp(-(map1.^2+map2.^2)/(2*sigma^2));

for i = 1+half_window_size:1:size(img_1,1)-half_window_size
    for j = size(img_1,2)-half_window_size:-1:1+half_window_size+d_max
    SAD_UNTILL_NOW = 1000;
    d_best = d_min;
    for d = d_min:1:d_max
        sample_1  = img_1(i-half_window_size:i+half_window_size,j-half_window_size:j+half_window_size);
        sample_2  = img_2(i-half_window_size:i+half_window_size,j-d-half_window_size:j-d+half_window_size);
        
        %SAD = sum(sum(abs(double(sample_1)-double(sample_2))));
        SAD1 = abs(double(sample_1)-double(sample_2));
        gauss_val = sum(sum(gauss*SAD1));

        if (gauss_val < SAD_UNTILL_NOW)
            SAD_UNTILL_NOW=gauss_val;
            d_best = d;
        end
    end
    disparity_map(i,j)=d_best;
   end
end

figure
imshow(disparity_map,[]);
colorbar
title('disparity map')
%% find depth map

% d = Mx*Tx*f / Z so Z = Mx*Tx*f/d

% suppose Tx = 100
Tx = 1
mx = 32
f = 31

depth_map = mx*f*Tx ./ disparity_map;

if (depth_map(:) == inf)
  depth_map(:) = nan;
end

figure
%imshow(depth_map,[0,0.1]); %thresholding for the second pair of images
imshow(depth_map,[1,80]); %thresholding for the second pair of images
colorbar
title ('depth map')