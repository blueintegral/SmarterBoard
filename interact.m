function interact()
%The main function to run and have magic happen.
close all;
%clear;
%clc;
%while(1)
   fprintf('Getting a new image!\n');
   %First fire off the python script to get an image from the camera
   %system('C:/Users/hunter/Desktop/mhacks/capture.py'); %disabled for demp
   %Creates an image called capture_raw.jpg that we now read in
   capture_raw = imread('capture_raw.jpg');
   imshow(capture_raw);
   %if (~capture_raw) fprintf('Failed to open raw image.');
   %Now we need to find the characters that we'll OCR later.
   %So we perform blob detection to find them. This assumes that the person
   %writing on the board has decent spacing between characters.
   capture_crop = imcrop(capture_raw, [160, 51 , 310, 260]);
   hold on;
   figure 
      imshow(capture_crop)

   capture_crop_gray = rgb2gray(capture_crop);
   
   %Now we have cropped out everything except the workspace on the board
   %Next we need to look for patches of dark pixels that we'll group 
   %into blobs
   
   %first divide it vertically. Look for the leftmost black pixel. That's
   %the first edge of the first number.
   %Then keep going until we hit a line with no black pixels, that's the
   %end of the first number. 
   %size(capture_crop)
   %301x305
   in_number = 0;
   in_number_y = 0;
   index = 0;
   cropped_images = {};
   for i = 1:310    
        if(capture_crop_gray(:, i, :) > 100)
            capture_crop(:, i, 1) = 255;  %red
            if (in_number)
                edge_back = i;
                in_number = 0;
                %Just finished traversing a number in the x direction, now
                %let's get the y coordinates
                %%%%%%%%%%%%%%%
                in_number_y = 0;
                for j = 1:260
                    if (capture_crop_gray(j, :, :) > 100)
                        capture_crop(j, :, 1) = 255; %red
                        if (in_number_y)
                            index = index + 1;  
                            %just finished traversing a number in the y
                            %direction
                            edge_bottom = j;
                            %store new image
                            cropped = capture_crop_gray(edge_top:edge_bottom, edge_front:edge_back);
                            %Now we need to resize them to be 16x16
                            size(cropped)
                            cropped = imresize(cropped, [16, NaN]);
                            [m, n] = size(cropped);
                            if(n < 16) %inflate to make 16x16
                                if(~mod(n, 2))
                                    %not divisible by 2, so we need pad with
                                    %1's unequally
                                    front_pad = (0.7*255).*ones(m, ((16-n)-1)/2);
                                    back_pad = (0.7*255).*ones(m, (((16-n)+1)/2)+1);
                                    %back_pad
                                    %front_pad
                                    cropped =  horzcat(front_pad, cropped);
                                    cropped = horzcat(cropped, back_pad);
                                else 
                                    %divisible by 2, pad with 1's equally
                                    front_pad = (0.7*255).*ones(m, (16-n)/2);
                                    back_pad = (0.7*255).*ones(m, ((16-n)/2)+1);
                                    cropped = horzcat(front_pad, cropped);
                                    cropped = horzcat(cropped, back_pad);
                                end
                            else
                                fprintf('bigger')
                                cropped = imresize(cropped, [16, NaN]);
                                [x y] = size(cropped);
                                step = y/16;
                                cropped = cropped(:, 1:step:end, :);
                                figure 
                                imshow(cropped)
                                
                            end
                            %cropped(:,:,:)
                            cropped_images{index} = cropped;
                            in_number_y = 0;
                        end
                    else
                        %we found the top of a character
                        if (~in_number_y)
                            edge_top = j;
                        end
                            in_number_y = 1;
                    end
                end
                %%%%%%%%%%%%%%%
            end
        else
            %we found the front of a character
            if (~in_number) 
                edge_front = i;
            end
                in_number = 1;
           
        end
       
   end
   imshow(capture_crop)
   fprintf('Number of objects found:\n');
   length(cropped_images)
   if (length(cropped_images) > 10)
        fprintf('OCR failed, attempting to fix...\n');
        interact()
   end
   for k = 1:length(cropped_images)
       %figure
       %size(cropped_images{k})
       %imshow(cropped_images{k})
       
       %needs a white number and black background
       level = graythresh(cropped_images{k});
       BW = im2bw(cropped_images{k},level);
       BW = imcomplement(BW);
       BW = circshift(BW, 1);
       imshow(BW)
       a3Test(BW)
       pause()
   end
%end

end