function F_best = F_ransac1(features1,features2)   
    
    f1 = [features1 ones(size(features1,1),1)];   % homogeneous coordinate
    f2 = [features2 ones(size(features2,1),1)];   % homogeneous coordinate        
    Max_error_accepted = 0.01;   
    iterations_ransan = 2000;
    max_C = 30;
    error_till_now =100000;
    
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
        
        t1 = [temp1 ones(size(temp1,1),1)];   % homogeneous coordinate
        t2 = [temp2 ones(size(temp2,1),1)];   % homogeneous coordinate
    
        if (m >= max_C)
            Ftemp = Ffinder(temp1,temp2);
            error_temp_model = sum(abs(sum((t2 .* (Ftemp * t1')'),2)));
            if ( error_temp_model <= error_till_now)
                 F_best = Ftemp;
                 error_till_now = error_temp_model;
            end
        end
    end
end