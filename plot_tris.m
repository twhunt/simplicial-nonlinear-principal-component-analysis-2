function tris_hndl = plot_tris(...
    tris, ...
    tri_coords_x, ...
    tri_coords_y, ...
    tri_coords_z, ...
    tris_hndl)

% colors =  mean(z(tris),2);
% colors = ((1+rd)/rd)*colors - rd;
if isempty(tris_hndl) || ~ishandle(tris_hndl);
    
    tris_hndl = trisurf(tris, tri_coords_x,tri_coords_y,tri_coords_z);
    
end


%set alpha channel so triangulated surface is translucent
%alpha(.8)
%set(tris_hndl, 'FaceColor', [.4 .4 .4]);
%colormap(bone);

num_tris = size(tris, 1);

%colormap(pink(256));
colormap(lines(256));
tmp_CData = get(tris_hndl, 'CData');

if isvector(tmp_CData)

    num_inds_to_gnrt = num_tris - numel(tmp_CData);
    tmp_rand_inds = round(256*rand(num_inds_to_gnrt, 1));
    tmp_CData = [tmp_CData(:); tmp_rand_inds];
    
else
   
   tmp_CData = round(256*rand(num_tris, 1));
    
end


set(tris_hndl, ...
    'CDataMapping', 'direct', ...
    'CData', tmp_CData, ...
    'FaceColor', 'flat',...
    'FaceAlpha', 1, ...
    'Faces', tris, ...
    'Vertices', [tri_coords_x(:), tri_coords_y(:) tri_coords_z(:)]);

