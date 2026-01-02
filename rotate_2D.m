function rotation_matrix = rotate_2D(theta) 
    % compute the 2D rotation matrix
    rotation_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
end