%% Obtain 2 cameras seeing the same scene

clear all;
close all;
clc;

PTS= load('pts3D.txt');
count_3d_points = length(PTS);
tmp = ones(count_3d_points,1);
PTS_homogenous = [PTS, tmp];


K1 = [8.0000000005798370e+002     2.5311355782387189e-008	  2.5600000001285036e+002;	
      0.0000000000000000e+000	  8.0000000007157689e+002	  2.5599999990332944e+002;	
      0.0000000000000000e+000	  0.0000000000000000e+000	  1.0000000000000000e+000];


C1 =  [1.2000000000891914e+003; 1.2000000000592427e+003; 1.2000000000697091e+003];  
  
theta1 = 180;
beta1 = 0;
alpha1 = 0;
Rz1 = [cosd(alpha1) -sind(alpha1) 0; sind(alpha1) cosd(alpha1) 0; 0 0 1];
Ry1 = [ cosd(theta1) 0 sind(theta1) ;0 1 0; -sind(theta1) 0 cosd(theta1)];
Rx1 = [ 1 0 0; 0 cosd(beta1) -sind(beta1); 0 sind(beta1) cosd(beta1)];

R1 = Rz1*Ry1*Rx1;

T1 = -R1*C1;  

Projection_Matrix1 = K1 * [R1, T1];

K2 = [8.0000000005798370e+002     2.5311355782387189e-008	  2.5600000001285036e+002;	
      0.0000000000000000e+000	  8.0000000007157689e+002	  2.5599999990332944e+002;	
      0.0000000000000000e+000	  0.0000000000000000e+000	  1.0000000000000000e+000];

theta2 = 45;
beta2 = 180;
alpha2 = 0;
Rz2 = [cosd(alpha2) -sind(alpha2) 0; sind(alpha2) cosd(alpha2) 0; 0 0 1];
Ry2 = [ cosd(theta2) 0 sind(theta2) ;0 1 0; -sind(theta2) 0 cosd(theta2)];
Rx2 = [ 1 0 0; 0 cosd(beta2) -sind(beta2); 0 sind(beta2) cosd(beta2)];

R2 = Rz2*Ry2*Rx2;

C2 =  [2.0000000000891914e+003; 1.0000000000592427e+003; 1.0000000000697091e+003];  
  
T2 = -R2*C2;  

Projection_Matrix2 = K2 * [R2, T2];

figure
X = PTS(:,1);
Y = PTS(:,2);
Z = PTS(:,3);

scatter3(X,Y,Z);
hold on

pose1 = rigid3d(R1,C1'); 
cam1=plotCamera('AbsolutePose', pose1, 'Opacity', 0, 'size', 50, 'Label', "First_Cam");
set(gca,'CameraUpVector',[0 -1 0]);

pose2 = rigid3d(R2,C2');
cam2=plotCamera('AbsolutePose', pose2, 'Opacity', 0, 'size', 50, 'Label', "Second_Cam");
set(gca,'CameraUpVector',[0 -1 0]);
hold off

%% Project 3D points with 2 different cameras

projected_points1 = zeros(count_3d_points, 3,1);
projected_points2 = zeros(count_3d_points, 3,1);
ProjectedX1 = zeros(count_3d_points, 2,1);
ProjectedY1 = zeros(count_3d_points, 2,1);

ProjectedX2 = zeros(count_3d_points, 2,1);
ProjectedY2 = zeros(count_3d_points, 2,1);

for i=1:1:count_3d_points
    projected_points1(i,:,:) =  Projection_Matrix1 * PTS_homogenous(i,:)';
    ProjectedX1(i) = projected_points1(i,1,:) ./ projected_points1(i,3,:);
    ProjectedY1(i) = projected_points1(i,2,:) ./ projected_points1(i,3,:);
end  

for i=1:1:count_3d_points
    projected_points2(i,:,:) =  Projection_Matrix2 * PTS_homogenous(i,:)';
    ProjectedX2(i) = projected_points2(i,1,:) ./ projected_points2(i,3,:);
    ProjectedY2(i) = projected_points2(i,2,:) ./ projected_points2(i,3,:);
end


figure
scatter(ProjectedX1(:),ProjectedY1(:));
hold on
scatter(ProjectedX2(:),ProjectedY2(:));
hold off



%% BackProject the 3D points using triangulation
coordinates=zeros([count_3d_points,3]);

for i=1:count_3d_points
    u = ProjectedX1(i);
    v = ProjectedY1(i);

    u_prime = ProjectedX2(i);
    v_prime = ProjectedY2(i);

    A(1, :) = v * Projection_Matrix1(3, :)' - Projection_Matrix1(2, :)';
    A(2, :) = u * Projection_Matrix1(3, :)' - Projection_Matrix1(1, :)' ;
    A(3, :) = v_prime * Projection_Matrix2(3, :)' - Projection_Matrix2(2, :)';
    A(4, :) = u_prime * Projection_Matrix2(3, :)' - Projection_Matrix2(1, :)';

    [U,S,V]=svd(A'*A)
    X=V(:, length(V))/V(end)
    
    % Add to coordinates matrix
    coordinates(i, 1)=X(1);
    coordinates(i, 2)=X(2);
    coordinates(i, 3)=X(3);
end
  
figure
scatter3(coordinates(:, 1), coordinates(:, 3), coordinates(:, 2));
hold on
pose1 = rigid3d(R1,C1'); 
cam1=plotCamera('AbsolutePose', pose1, 'Opacity', 0, 'size', 50, 'Label', "First_Cam");
set(gca,'CameraUpVector',[0 -1 0]);

pose2 = rigid3d(R2,C2');
cam2=plotCamera('AbsolutePose', pose2, 'Opacity', 0, 'size', 50, 'Label', "Second_Cam");
hold off
