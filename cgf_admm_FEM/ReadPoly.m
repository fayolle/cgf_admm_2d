function segments = ReadPoly(poly_name)
fid = fopen(poly_name, 'r');

% header
[line, ~] = fscanf(fid, '%d', 4);
% format is: # of vert., dim (2), # of attr., # of boundary markers
num_vertices = line(1);
% This .poly file is produced by Triangle and therefore the num of vert
% is always going to be 0 (i.e. the vertices are given in a separate file)
assert(num_vertices == 0);

[line, ~] = fscanf(fid, '%d', 2);
% format is: # of segments, # of boundary markers
num_segments = line(1);
num_boundary_markers = line(2);


% read segments line by line
segments = zeros(num_segments, 2);
one_start = true;
for i = 1:num_segments
    if num_boundary_markers
        [line, ~] = fscanf(fid, '%d', 4);
        idx = line(1);
        v1 = line(2);
        v2 = line(3);
    else
        [line, ~] = fscanf(fid, '%d', 3);
        idx = line(1);
        v1 = line(2);
        v2 = line(3);
    end


    % are arrays starting at 0 or 1?
    if (i == 1) && (idx == 0)
        one_start = false;
    end

    segments(i, 1) = v1;
    segments(i, 2) = v2;
end

% the last lines are:
% - the number of holes and the coordinates of a point in each hole
% -  some comments
% ignore them

% matlab uses array indices starting at 1
if (one_start == false)
    segments(:) = segments(:) + 1;
end


fclose(fid);
end
