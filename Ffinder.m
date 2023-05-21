function finalF = Ffinder(features1, features2)

    % data normalization 
    mean_x1 = mean(features1(:,1)); 
    mean_y1 = mean(features1(:,2));
    mean_x2 = mean(features2(:,1));
    mean_y2 = mean(features2(:,2));
    
    std_x1 = std(features1(:,1));
    std_x2 = std(features2(:,1));
    std_y1 = std(features1(:,2));
    std_y2 = std(features2(:,2));
    
    
    M1 = [1/(std_x1) 0 0; 0 1/(std_y1) 0; 0 0 1] * [1 0 -mean_x1; 0 1 -mean_y1; 0 0 1];
    M2 = [1/(std_x2) 0 0; 0 1/(std_y2) 0; 0 0 1] * [1 0 -mean_x2; 0 1 -mean_y2; 0 0 1];

    N = size(features1,1);
    % apply normalization
    for i = 1:1:N
        a = M1 * [features1(i,1); features1(i,2); 1];
        x1_norm(i) = a(1); y1_norm(i) = a(2);
        b = M2 * [features2(i,1); features2(i,2); 1];
        x2_norm(i) = b(1); y2_norm(i) = b(2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fundamantal matrix -> called A at this step (Ax = 0)
    if N == 8 
        A = zeros(9,9);
    end
    if N ~= 8
        A = zeros(N,9);
    end

    for i=1:1:N
       A(i,:) = [(x1_norm(i))*(x2_norm(i)) (y1_norm(i))*(x2_norm(i)) (x2_norm(i)) (x1_norm(i))*(y2_norm(i)) ...
               (y1_norm(i))*(y2_norm(i)) (y2_norm(i)) (x1_norm(i)) (y1_norm(i)) 1];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply SVD (on A)
    [U,S,V] = svd(A);
    
    % find F - the eigenvector that corresponds to the smallest eigenvalue
    F_norm_rank3 = [V(1,9),V(2,9),V(3,9);...
        V(4,9), V(5,9),V(6,9);...
        V(7,9), V(8,9),V(9,9)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rank(F_rank3) = 3 , so we need to make it have rank 2
    [UF,SF,VF] = svd(F_norm_rank3);
    SF(3,3) = 0 ;
    F_norm_rank2 = UF*SF* transpose(VF); % This F has rank 2 
    % denormalize F
    finalF  = transpose(M2) * F_norm_rank2 * M1;
    finalF = finalF/norm(finalF);
end