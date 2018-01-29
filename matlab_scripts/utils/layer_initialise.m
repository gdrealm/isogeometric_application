function layer_initialise(fid, str)

fprintf(fid, '# geometry layers data file for isogeometric application\n');
%fprintf(fid, '# (c) %d Hoang Giang Bui, Ruhr-University Bochum\n',year(now));
fprintf(fid, '# (c) %d Hoang Giang Bui, Ruhr-University Bochum\n',2018);
c = clock;
fprintf(fid, '# This file is created at %s\n',datestr(datenum(c(1),c(2),c(3),c(4),c(5),c(6))));
fprintf(fid, '%s\n', str);
fprintf(fid, 'class LayerProvider():\n');
fprintf(fid, '\tdef __init__(self):\n');
fprintf(fid, '\t\tself.layer_list = []\n');
% fprintf(fid, '\t\t#self.layer_list provides the list of layer in the model\n');
fprintf(fid, '\t\tself.layer_attributes = {}\n');
fprintf(fid, '\t\tself.layer_nodes_sets = {}\n');
fprintf(fid, '\t\tself.layer_entities_sets = {}\n');
fprintf(fid, '\t\tself.layer_entity_info_sets = {}\n');
fprintf(fid, '\t\tself.layer_boundary_marker = {}\n');
fprintf(fid, '\t\tself.control_points_row_u = {}\n');
fprintf(fid, '\t\tself.control_points_row_v = {}\n');
fprintf(fid, '\n');

end
