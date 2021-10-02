clear all;
close all;
clc;

PTS1 = load('pts2D_1.txt'); % 300 x 2 matrix
PTS2 = load('pts2D_2.txt'); % 300 x 2 matrix

tmp = ones(length(PTS1),1);

Homogenous_PTS1 = [PTS1, tmp];
Homogenous_PTS2 = [PTS2, tmp]; % 300 x 3 matrix (z axis = 1 for each)
 
matches = [PTS1, PTS2]; % 300 x 4 matrix. Firs two row for x-y values of first image, second two row for x-y values of second image
 
%% Create unnormalized Fundamental Matrix

A=ones(size(matches,1),9);

A(:,1)=matches(:,1).*matches(:,3);
A(:,2)=matches(:,1).*matches(:,4);
A(:,3)=matches(:,1);
A(:,4)=matches(:,2).*matches(:,3);
A(:,5)=matches(:,2).*matches(:,4);
A(:,6)=matches(:,2);
A(:,7)=matches(:,3);
A(:,8)=matches(:,4);

%first SVD
[~,S,V]=svd(A);
%find min singular value
[~,I]=min(diag(S));
f=V(:,I);

disp("unnormalized fundamental matrix")
F_matrix = reshape(f,3,3)'

 %% Enforcing rank 2 constraint 
 
[U1,S1,V1]=svd(F_matrix);
 
S1

[~,I1]=min(diag(S1))
 
S1(I1,I1)=0; % but its already 0 so actually nothing changes. so this constraint is already satisfied in our case
 
disp("unnormalized but rank 2 constraint applied fundamental matrix")
F_matrix = U1*S1*V1'
 
 
 %% Calculating and Drawing Epipolar Lines with Fundamental Matrix
 
 Epipolar_Lines1 = ones([length(PTS1),3]);
 Epipolar_Lines2 = ones([length(PTS2),3]);
 
 x1 = ones([length(PTS1),1]);
 x2 = ones([length(PTS1),512]);
  
 for i=1:length(PTS1)
     Epipolar_Lines1(i,:) = F_matrix * Homogenous_PTS2(i,:)';
     Epipolar_Lines2(i,:) = F_matrix' * Homogenous_PTS1(i,:)';
 end
 
 figure
 title("epipolar lines of 1. image with unnormalized fundamental matrix")
 % Compute the intersection points of the lines and the image border.
 scatter(PTS1(:,1), PTS1(:,2));
 hold on
 pts = lineToBorderPoints(Epipolar_Lines1, [512,512]);
 % Show the epipolar lines in the first image
 line(pts(:, [1,3])', pts(:, [2,4])');
 
 figure
 title("epipolar lines of 2. image with unnormalized fundamental matrix")
 scatter(PTS2(:,1), PTS2(:,2));
 hold on
 pts = lineToBorderPoints(Epipolar_Lines2, [512,512]);
 line(pts(:, [1,3])', pts(:, [2,4])');
 
 % compute and draw just for 1 point. To see what is going on in fact
 figure
 title("First epipolar line of 1. image with unnormalized fundamental matrix")
 scatter(PTS1(1,1), PTS1(1,2));
 hold on
 pts = lineToBorderPoints(Epipolar_Lines1, [512,512]);
 line(pts(1, [1,3])', pts(1, [2,4])');
 hold off
 
 figure
 title("First epipolar line of 2. image with unnormalized fundamental matrix")
 scatter(PTS2(1,1), PTS2(1,2));
 hold on
 pts = lineToBorderPoints(Epipolar_Lines2, [512,512]);
 line(pts(1, [1,3])', pts(1, [2,4])');
 hold on


error1 = zeros(length(PTS1),1);
error2 = zeros(length(PTS1),1);
 
% computing average error
for i=1:length(PTS1)
    error1(i,:) = Homogenous_PTS1(i,:) * Epipolar_Lines1(i,:)';
    error2(i,:) = Homogenous_PTS2(i,:) * Epipolar_Lines2(i,:)';
end

average_error = (sum(error1,'all') + sum(error2,'all')) / 2 

 %% Intersections of epipolar lines (epipoles)

epipole1 = null(F_matrix')
epipole2 = null(F_matrix)



%%  Extraction the camera projection matrices from the Fundamental matrix

epipole2_skew_symmetric = [0 -epipole2(3) epipole2(2)
                           epipole2(3) 0  -epipole2(1) 
                          -epipole2(2) epipole2(1)  0];


P2 = [epipole2_skew_symmetric*F_matrix, epipole2]
P1 = [eye(3), [0;0;0]]
