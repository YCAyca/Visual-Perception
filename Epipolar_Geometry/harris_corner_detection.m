% reference : https://medium.com/data-breach/introduction-to-harris-corner-detector-32a88850b3f6
clear all; 
close all; 
clc;


% Take the grayscale of the original image

im=imread('corner2.png');
figure
imshow(im)
im = rgb2gray(im); % its already grayscale but better for work different images

% Apply a Gaussian filter to smooth out any noise - decided to not to use on my images
% im = imgaussfilt (im, 4);
% figure
% imshow(im)

% Apply Sobel operator to find the x and y gradient values for every pixel in the grayscale image
dx=[1 0 -1; 2 0 -2; 1 0 -1];    
dy=[1 2 1; 0 0 0; -1 -2 -1];     
 
Ix=conv2(im, dx, 'same');
Iy=conv2(im, dy, 'same');

% For each pixel p in the grayscale image, 
% consider a 3×3 window around it and compute the corner strength function. Call this its Harris value

Ix2 = Ix .^ 2;
Iy2 = Iy .^ 2;
Ixy = Ix .* Iy;

kernel_size = 3;

% Compute the sums of the products of derivatives at each pixel according
% to the kernel size

sum_all = ones([kernel_size,kernel_size]);

Sx2 = conv2(Ix2, sum_all, 'same');
Sy2 = conv2(Iy2, sum_all, 'same');
Sxy = conv2(Ixy, sum_all, 'same');

kernel_size = 3;
k = 0.04;
Threshold = 2.5*10e+09;

[numOfRows, numOfColumns] = size(im);

R_avg = 0;
for x=1:numOfRows
   for y=1:numOfColumns
       % 4) Define at each pixel(x, y) the matrix M
       M = [Sx2(x, y) Sxy(x, y); Sxy(x, y) Sy2(x, y)];
       
       % 5) Compute the response of the detector at each pixel
       R = det(M) - k * (trace(M) ^ 2);
       
       % 6) Threshold on value of R
       if (R > Threshold)
          im = insertMarker(im,[y,x]);
          
       %take sum of R for average after
       R_avg = R_avg + R;
       end
   end
end
  

R_avg = R_avg / (numOfRows*numOfColumns);

disp(sprintf("average of R is : %d", R_avg));
disp(sprintf("So its better to try threshold between: %d", R_avg * 4));
disp(sprintf("and %d" , R_avg * 5));

figure
imshow(im)

im=imread('corner2.png');
im = rgb2gray(im);
%im = imgaussfilt (im, 4);
C = corner(im);
figure
imshow(im)
hold on
plot(C(:,1),C(:,2),'r*');

       
