function ring = compute_vertex_face_ring(face)
% compute_vertex_face_ring - compute the faces adjacent to each vertex

nfaces = size(face,1);
nverts = max(face(:));

ring{nverts} = [];

for i=1:nfaces
    for k=1:3
        ring{face(i,k)}(end+1) = i;
    end
end
end
