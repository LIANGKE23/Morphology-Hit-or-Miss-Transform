% Purpose::
%        to explore mathematical morphology and use the 
%        hit-or-miss transform to detect the smallest and 
%        largest circles in an image.
% 
% Processing Flow::
%        1. Pre-process the image- create binary image and filter to remove
%           salt and pepper noise
%        2. Design the structuring element
%        3. Detect the smallest and largest circles
% Authors::
%       Lindsey Schwartz, Ke Liang, Xilun Liu
% Date::
%       01/22/2019

%% Load image and pre-process
% read in the image and get the dimensions
im = imread('F:\PSU第一学期相关书籍及其工作\Digtal Image Processing II\Project 1\Project1\Project1\RandomDisks-P10.jpg');
[num_rows, num_cols, num_bands] = size(im);
figure, imshow(im), title('Original image')

% convert to black and white
if num_bands > 1
    im = rgb2gray(im);
end
bw_im = zeros(size(im(:,:,1)));
bw_thresh = 128;
for row = 1:num_rows
    for col = 1:num_cols
        if im(row,col,1) >= bw_thresh
            bw_im(row,col) = 1;
        end
    end
end

figure, imshow(bw_im), title('Image converted to black and white')

% remove salt and pepper noise
filt_im = removeNoise(bw_im, num_rows, num_cols);

figure, imshow(filt_im), title('Salt and pepper noise removed')


%% Hit-or-Miss Transform
% Pick the structuring element (by any means necessary)
[kernel_sm_A, kernel_sm_B, kernel_lg_A, kernel_lg_B, min_groups, max_groups] = getStructuringElement(filt_im);

figure,
h1 = subplot(2,2,1);
imshow(kernel_sm_A)
title('Small structuring element A')
h2 = subplot(2,2,2);
imshow(kernel_sm_B)
title('Small structuring element B')
h3 = subplot(2,2,3);
imshow(kernel_lg_A)
title('Large structuring element A')
h4 = subplot(2,2,4);
imshow(kernel_lg_B)
title('Large structuring element B')
linkaxes([h4,h3,h2,h1])

% move the structuring element through the image to detect hits
hits_small = zeros(size(filt_im));
r = ((length(kernel_sm_A(:,1))-1)/2);
for i = r + 1:num_rows - r - 1
    for j = r + 1:num_cols - r - 1
        nhood = filt_im(i-r:i+r,j-r:j+r);
        response = sum(sum(kernel_sm_A + nhood));
        if response == 0
            hits_small(i,j) = 1;
        end
    end
end

hits_large = zeros(size(filt_im));
r = ((length(kernel_lg_A(:,1))-1)/2);
for i = r + 1:num_rows - r - 1
    for j = r + 1:num_cols - r - 1
        nhood = filt_im(i-r:i+r,j-r:j+r);
        response = sum(sum(kernel_lg_A + nhood));
        if response == 0
            hits_large(i,j) = 1;
        end
    end
end

% move the structuring element through the image to detect misses
r = ((length(kernel_sm_B(:,1))-1)/2);
misses_small = zeros(size(filt_im));
complement_im = ~filt_im;
for i = r + 1:num_rows - r - 1
    for j = r + 1:num_cols - r - 1
        nhood = complement_im(i-r:i+r,j-r:j+r);
        nhood_sum = sum(nhood(:));
        response = kernel_sm_B.*nhood;
        response_sum = sum(response(:));
        if response_sum >= nhood_sum
            misses_small(i,j) = 1;
        end
    end
end

r = ((length(kernel_lg_B(:,1))-1)/2);
misses_large = zeros(size(filt_im));
complement_im = ~filt_im;
for i = r + 1:num_rows - r - 1
    for j = r + 1:num_cols - r - 1
        nhood = complement_im(i-r:i+r,j-r:j+r);
        nhood_sum = sum(nhood(:));
        response = kernel_lg_B.*nhood;
        response_sum = sum(response(:));
        if response_sum >= nhood_sum
            misses_large(i,j) = 1;
        end
    end
end

% smallest and largest circles detected
small = hits_small.*misses_small;
large = hits_large.*misses_large;

% Plot the smallest and largest circles on the original image
detected_circles = ones(size(im));
r = mean([min_groups(1) max_groups(1)]);
figure, imshow(im)
hold on
for i = 1:num_rows
    for j = 1:num_cols
        if small(i,j) == 1
            ang=0:0.01:2*pi; 
            xp=r*cos(ang);
            yp=r*sin(ang);
            plot(j+xp,i+yp, 'Color', 'r','LineWidth',1);
            for k = 1:num_rows
                for l = 1:num_cols
                    dist = sqrt((k-i)^2+(l-j)^2);
                    if dist <= r
                        detected_circles(k,l) = 0;
                    end
                end
            end
        end
    end
end

r = mean([min_groups(end) max_groups(end)]);
for i = 1:num_rows
    for j = 1:num_cols
        if large(i,j) == 1
            ang=0:0.01:2*pi; 
            xp=r*cos(ang);
            yp=r*sin(ang);
            plot(j+xp,i+yp, 'Color', 'r','LineWidth',1);
            for k = 1:num_rows
                for l = 1:num_cols
                    dist = sqrt((k-i)^2+(l-j)^2);
                    if dist <= r
                        detected_circles(k,l) = 0;
                    end
                end
            end
        end
    end
end
hold off

figure, imshow(detected_circles), title('Smallest and largest detected circles')
