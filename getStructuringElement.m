function [kernel_sm_A, kernel_sm_B, kernel_lg_A, kernel_lg_B, min_groups, max_groups] = getStructuringElement(filt_im)
% Purpose::
%       To automatically get the structuring elements to find the smallest
%       and largest circles in the image
% Inputs::
%       filt_im - a binary image with noise removed
% Outputs::
%       kernel_sm_A - smallest structure to match
%       kernel_sm_B - smallest structure that doesn't match
%       kernel_lg_A - largest structure to match
%       kernel_lg_B - largest structure that doesn't match


%% Get the structuring elements
% determine the size of the circles
stats = regionprops('table',~filt_im,'Centroid','MajorAxisLength','MinorAxisLength');
diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii = diameters/2;

% initialize the classes
group1 = [];
group2 = [];
group3 = [];
group4 = [];
group5 = [];

% sort the circles into 5 classes(sizes)
idx = kmeans(radii(2:end),5);
for i = 1:numel(idx)
    if idx(i) == 1
        group1 = [group1 radii(i+1)];
    end
    if idx(i) == 2
        group2 = [group2 radii(i+1)];
    end
    if idx(i) == 3
        group3 = [group3 radii(i+1)];
    end
    if idx(i) == 4
        group4 = [group4 radii(i+1)];
    end
    if idx(i) == 5
        group5 = [group5 radii(i+1)];
    end
end

% determine the size of the smallest and largest circles
min_groups = sort([min(group1) min(group2) min(group3) min(group4) min(group5)]);
max_groups = sort([max(group1) max(group2) max(group3) max(group4) max(group5)]);

r_small_min = floor(min_groups(1));
r_small_max = floor(max_groups(1));

r_large_min = ceil(min_groups(end));
r_large_max = ceil(max_groups(end));

% create the structuring element based on the smallest and largest circles
% the structuring element is a square
% largest square that can fit inside a circle is d/sqrt(2)
kernel_sm_A = zeros([floor(r_small_min*2/sqrt(2))+1 floor(r_small_min*2/sqrt(2))+1]);
kernel_sm_B = ones([ceil(r_small_max*2)+1 ceil(r_small_max*2)+1]);
kernel_sm_B = padarray(kernel_sm_B,[5 5],0,'both');

kernel_lg_A = zeros([floor(r_large_min*2/sqrt(2)) floor(r_large_min*2/sqrt(2))]);
kernel_lg_B = ones([ceil(r_large_max*2)+1 ceil(r_large_max*2)+1]);
kernel_lg_B = padarray(kernel_lg_B,[5 5],0,'both');
