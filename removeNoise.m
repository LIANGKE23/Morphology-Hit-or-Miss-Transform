function [filt_im] = removeNoise(bw_im, num_rows, num_cols)
% Purpose::
%       To remove salt and pepper noise from an image
% Inputs::
%       bw_im - a binary image
%       num_rows - number of rows in bw_im
%       num_cols - number of columns in bw_im
% Outputs::
%       filt_im - the image with salt and pepper noise removed


%%
filt_im = zeros(size(bw_im));
for i = 2:num_rows - 1
    for j = 2:num_cols - 1
        nhood = [bw_im(i-1,j-1), bw_im(i-1,j), bw_im(i-1,j+1); bw_im(i,j-1), bw_im(i,j), bw_im(i,j+1); bw_im(i+1,j-1), bw_im(i+1,j), bw_im(i+1,j+1)];
        filt_im(i,j) = median(nhood(:));
    end
end


%%
% % open - 1 pixel
% open_im = bw_im;
% kernel = [1 0 1];
% for row = 1:num_rows
%     for col = 2:num_cols - 1
%         nhood = [bw_im(row,col-1), bw_im(row,col), bw_im(row,col+1)];
%         union = max(nhood,kernel);
%         if nhood == union
%             open_im(row,col) = 1;
%         end
%     end
% end
% 
% % close - 1 pixel
% close_im = open_im;
% kernel = [0 1 0];
% for row = 1:num_rows
%     for col = 2:num_cols - 1
%         nhood = [open_im(row,col-1), open_im(row,col), open_im(row,col+1)];
%         union = min(nhood,kernel);
%         if nhood == union
%             close_im(row,col) = 0;
%         end
%     end
% end
% 
% % open - 2 pixels
% kernel = [1 0 0 1];
% for row = 1:num_rows
%     for col = 2:num_cols - 2
%         nhood = [close_im(row,col-1), close_im(row,col), close_im(row,col+1) close_im(row,col+2)];
%         union = max(nhood,kernel);
%         if nhood == union
%             close_im(row,col) = 1;
%         end
%     end
% end
% 
% % close - 2 pixels
% kernel = [0 1 1 0];
% for row = 1:num_rows
%     for col = 2:num_cols - 2
%         nhood = [close_im(row,col-1), close_im(row,col), close_im(row,col+1) close_im(row,col+2)];
%         union = min(nhood,kernel);
%         if nhood == union
%             close_im(row,col) = 0;
%         end
%     end
% end
% filt_im = close_im;