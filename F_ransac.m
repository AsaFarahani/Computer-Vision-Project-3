function F_best = F_ransac(features1,features2)   

    f1 = [features1 ones(size(features1,1),1)];   % homogeneous coordinate
    f2 = [features2 ones(size(features2,1),1)];   % homogeneous coordinate
    j = 0;
    Max_number_inliers = 0;         
    Max_error_accepted = 0.05;   
    iterations_ransan = 10000;

    for i = 1:1:iterations_ransan    
        % randomly select 8 pairs so that you would be able to find F 
        sampled_8_indices = randsample(size(features1,1),8);
        features1_8points = features1(sampled_8_indices,:);
        features2_8points = features2(sampled_8_indices,:);
        
        % use a function called (Ffinder) to find the normF - it is rank 2 
        % exact method
        F = Ffinder(features1_8points, features2_8points);
    
        error = sum((f2 .* (F * f1')'),2);
        Inliers = size(find(abs(error) <= Max_error_accepted) , 1);
    
        m = 0;
        for g = 1:1:size(features1,1)  
            if (abs(error(g,:)) <= Max_error_accepted)
                m = m + 1;
                temp1(m,:) = features1(g,:);
                temp2(m,:) = features2(g,:);
            end
        end
    
        if (Inliers > Max_number_inliers)
           j = j+1;
           F_best = Ffinder(temp1,temp2);    
           Max_number_inliers = Inliers;
        end 
    end
end