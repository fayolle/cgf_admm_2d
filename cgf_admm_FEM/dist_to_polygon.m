function d = dist_to_polygon(nodes, segments)
% unsigned distance for each point in nodes to the mesh
% defined by (nodes, poly)

% ignore the last columns if any
nodes = nodes(:,1:2);

nnodes = size(nodes, 1);
d = zeros(nnodes, 1);
for i=1:nnodes
    p = nodes(i,:);
    d(i) = dist_to_polygon_from_point(p, nodes, segments);
end

end


function d = dist_to_polygon_from_point(p, nodes, segments)
% unsigned distance from point p to the polygon described by (nodes, poly)

npoly = size(segments, 1);
dist = Inf;

for i=1:npoly
    v1_idx = segments(i,1);
    v2_idx = segments(i,2);
    v1 = nodes(v1_idx,:);
    v2 = nodes(v2_idx,:);

    cdist = dist_to_segment(p, v1, v2);
    dist = min(dist, cdist);
end

d = dist;

end


function dist = dist_to_segment(p, v1, v2)
param = dotp(p-v1, v2-v1);
param = param / squared_norm(v2-v1);
if param < 0
    dist = euclidean_dist(v1, p);
elseif param > 1
    dist = euclidean_dist(v2, p);
else
    dist = euclidean_dist(p, v1+param*(v2-v1));
end
end



function sqn = squared_norm(v)
% return norm(v)^2 where norm(v) is Sqrt[x^2 + y^2], v = [x,y]
sqn = v(1)^2 + v(2)^2;
end


function dist = euclidean_dist(p1, p2)
% return the Euclidean distance between the 2D points p1 and p2
dist = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
end


function s = dotp(v1, v2)
s = v1(1)*v2(1) + v1(2)*v2(2);
end


