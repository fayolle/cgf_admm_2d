function ele = ReadElements(ele_name)
fid = fopen(ele_name, "r");
% read the header
[line, ~] = fscanf(fid, '%d', 3);
num_ele = line(1);

ele = zeros(num_ele, 3);
one_start = true;
for i = 1:num_ele
    [line, ~] = fscanf(fid, '%d', 4);
    idx = line(1);
    v1 = line(2);
    v2 = line(3);
    v3 = line(4);

    % are arrays starting at 0 or 1?
    if (i == 1) && (idx == 0)
        one_start = false;
    end

    ele(i, 1) = v1;
    ele(i, 2) = v2;
    ele(i, 3) = v3;
end

% the last lines are some comments; ignore them

% matlab uses array indices starting at 1
if (one_start == false)
    ele(:) = ele(:) + 1;
end

fclose(fid);
end
