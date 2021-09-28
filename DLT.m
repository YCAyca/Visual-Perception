% Apply DLT using the following functions 
% World_PTS : 3D world points, GT_PTS: 2D Projected image points of these 3D world points, 
% Calibration_Points_Count : the quantity of points we use for realizing calibration

function DLT(World_PTS, GT_PTS, Calibration_Points_Count)
   A = Create_Matrix_A(Calibration_Points_Count, World_PTS, GT_PTS)
   P = Obtain_Projection_Matrix(A)
   [Projected_X, Projected_Y] = Projection(P, World_PTS, GT_PTS);
end

function [A] = Create_Matrix_A(total_number_of_points, World_PTS, GT_PTS)
    World_PTS = World_PTS(1:total_number_of_points,:);
    X = World_PTS(:, 1);
    Y = World_PTS(:, 2);
    Z = World_PTS(:, 3);
    GT_PTS = GT_PTS(1:total_number_of_points,:); % Ground truth points
    x = GT_PTS(:, 1);
    y = GT_PTS(:, 2);

    % Creating matrix A with giving points 
    for i=1:total_number_of_points
        A(2*i-1,:) = [-X(i) -Y(i) -Z(i) -1 0 0 0 0 x(i)*X(i) x(i)*Y(i) x(i)*Z(i) x(i)] ;
        A(2*i,:) = [0 0 0 0 -X(i) -Y(i) -Z(i) -1 y(i)*X(i) y(i)*Y(i) y(i)*Z(i) y(i)] ;
    end
end

function [P] = Obtain_Projection_Matrix(A)
    % Compute P from A using SVD
    [U,S,V] = svd(A'*A);
    
    unknown_parameters = V(:,end);

    P_transpoze = reshape(unknown_parameters, [4,3]);
    P = P_transpoze';
end

% Function for testing the result! 
% This time we use the projection matrix we found to calculate 2D image points from the same 3D world points we use during calibration
function [Projected_X, Projected_Y] = Projection(P, World_PTS, GT_PTS)
    total_3d_points = length(GT_PTS); % ground truth 2D image plane points
    tmp = ones(total_3d_points,1);
    projected_points = zeros(total_3d_points, 3,1); % 300 x 3 x 1 matrix
    homogenous_points = [World_PTS, tmp]; % world 3D points converted homogenous points

    for i=1:1:total_3d_points
        projected_points(i,:,:) =  P * homogenous_points(i,:)'; 
    end   

    figure

    scatter(GT_PTS(:,1),GT_PTS(:,2),'filled');
    hold on

    figure
    Projected_X = projected_points(:,1,:) ./ projected_points(:,3,:);
    Projected_Y = projected_points(:,2,:) ./ projected_points(:,3,:);

    scatter(Projected_X,Projected_Y);
    hold off
end
