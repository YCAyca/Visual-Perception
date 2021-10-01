clear all;
close all;
clc;

%We have simulated a regular checkerboard pattern, comprising 10x10 points, and a camera observing it. 
%The checkerboard is moved around the camera and 10 images of it are acquired in the process. 

%PTS_world = load('ptsXY.txt'); % Cartesian coordinates of the checkerboard points,
%Pixel coordinates of the 2D points in the 10 images, given in pts2D_i.txt, i = 1, 2, ..., 10.

%PTS_image_plane(:,:,1); % 100 2D points coming from first image 

%% Obtaining 10 homography matrices H for 10 different images

number_of_points = 100;
number_of_images = 10;

PTS_world = load('ptsXY.txt');

PTS_image_plane = zeros([2,number_of_points,number_of_images]);

for i=1:1:number_of_images
   file_name = "pts2D_" + int2str(i) + ".txt";
   PTS_image_plane(:,:,i) = load(file_name);
end

A = zeros(number_of_points,9);
H = zeros([3,3,10]); % 10 times 3x3 homography matrix

for k=1:1:number_of_images % 10 different images, 10 different homography matrixes!
    for i=1:1:number_of_points
        x = PTS_image_plane(1,i,k);
        y = PTS_image_plane(2,i,k);
        X = PTS_world(i,1); 
        Y = PTS_world(i,2); 
        A(2*i-1,:) = [X Y 1 0 0 0 -x*X -x*Y -x];
        A(2*i,:) = [0 0 0 X Y 1 -y*X -y*Y -y];
    end
    [U,D,V] = svd(A);
    % Extract homography
    H(:,:,k) = reshape(V(:,9),3,3)';
end


V = zeros([2*10,6]);
for i=1:1:number_of_images
    v12_trans  = [ H(1,1,i)*H(1,2,i);
                   H(1,1,i)*H(2,2,i) + H(2,1,i)*H(1,2,i)
                   H(2,1,i)*H(2,2,i)
                   H(3,1,i)*H(1,2,i) + H(1,1,i)*H(3,2,i)
                   H(3,1,i)*H(2,2,i) + H(2,1,i)*H(3,2,i)
                   H(3,1,i)*H(3,2,i)]';  

    v11  = [H(1,1,i)*H(1,1,i)
            H(1,1,i)*H(2,1,i) + H(2,1,i)*H(1,1,i)
            H(2,1,i)*H(2,1,i)
            H(3,1,i)*H(1,1,i) + H(1,1,i)*H(3,1,i)
            H(3,1,i)*H(2,1,i) + H(2,1,i)*H(3,1,i)
            H(3,1,i)*H(3,1,i)];


    v22  = [H(1,2,i)*H(1,2,i)
            H(1,2,i)*H(2,2,i) + H(2,2,i)*H(1,2,i)
            H(2,2,i)*H(2,2,i)
            H(3,2,i)*H(1,2,i) + H(1,2,i)*H(3,2,i)
            H(3,2,i)*H(2,2,i) + H(2,2,i)*H(3,2,i)
            H(3,2,i)*H(3,2,i)];
        
    
     V(2*i-1,:) =  v12_trans;
     V(2*i,:) = (v11-v22)';
end

size(V)

[~,~,B] = svd(V);
    
b = B(:,end); % b is 6D vector


%% Camera intrinsic parameters (K)

x0 = (b(2)*b(4)-b(1)*b(5)) / (b(1)*b(3)-b(2)^2);
lambda = b(6) - (b(4)^2 + x0*(b(2)*b(4)-b(1)*b(5))) / b(1);


fx = sqrt(lambda / b(1));
fy = sqrt(lambda*b(1) / (b(1)*b(3)-b(2)^2));

gamma = -b(2)*fx^2*fy/lambda;
y0 = gamma*x0 / fx - b(4)*fx^2 / lambda;


K = [fx 0 x0; 0 fy y0; 0 0 1]
% check with ground truth K matrix
K_groundtruth = load('K.txt')

%% Camera extrinsic parameters (R,t) -just the last H matrix used, not calculated for 10 different H -

r1 = lambda * (K\H(:,1,1)); % \ = inv(K) *
r2 = lambda * (K\H(:,2,1));
r3 = cross(r1, r2);

R = [r1 r2 r3]
t = lambda * (K\H(:,3,1))
