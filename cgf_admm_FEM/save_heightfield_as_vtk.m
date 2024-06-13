% Very preliminary and untested code
% Save the poisson_dist field to a file using the (legacy) VTK format


function save_heightfield_as_vtk(nodes, elements, field, output_filename, title)
%
% writes the result of a FEM computation to a legacy vtk file
%

assert(size(nodes,1) ~= 0);
assert(size(elements,1) ~= 0);
assert(size(field,1) ~= 0);
% other arguments are optional


[number_nodes, ~] = size(nodes);
[number_elements, element_order] = size(elements);
% support only linear element (in 2D: triangles)
assert(element_order == 3);


if (isempty(output_filename))
    output_filename = 'solution.vtk';
end

if (isempty(title))
    title = 'solution';
end


xyz = zeros(3, number_nodes);
xyz(1:2, :) = nodes(:,1:2)';
xyz(3, :) = field(:);


% matlab indices are starting at 1;
% in vtk they start from 0
element_node = zeros(element_order, number_elements);
element_node(1:element_order, :) = elements(:, 1:element_order)' - 1;


fid = fopen(output_filename, 'w');

vtk_write_field(fid, title, number_nodes, number_elements, element_order, xyz, element_node, field');

fclose(fid);

end


function vtk_write_field(fid, title, number_nodes, number_elements, element_order, xyz, element_node, field)

fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, '%s\n', title);
fprintf(fid, 'ASCII\n');
fprintf(fid, '\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid, 'POINTS %d double\n', number_nodes);


% output nodes
for i = 1:number_nodes
    fprintf(fid, '%f %f %f\n', xyz(1:3, i));
end


cell_size = number_elements * (element_order + 1);
fprintf(fid, '\n');
fprintf(fid, 'CELLS %d %d\n', number_elements, cell_size);
for i = 1:number_elements
    fprintf(fid, ' %d', element_order);
    for j = 1:element_order
        fprintf(fid, ' %d', element_node(j, i));
    end
    fprintf(fid, '\n');
end


fprintf(fid, '\n');
fprintf(fid, 'CELL_TYPES %d\n', number_elements);

assert(element_order==3);

% linear elements (triangle in 2D) is cell type 5
for i = 1:number_elements
    fprintf(fid, '5\n');
end


% save the field either defined at each node or at each element
field = full(field);

fprintf(fid, '\n');
fprintf(fid, 'POINT_DATA %d\n', number_nodes);
fprintf(fid, 'SCALARS field double\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
for i = 1:number_nodes
    fprintf(fid, ' %f\n', field(1, i));
end


fprintf(fid, '\n');
end

