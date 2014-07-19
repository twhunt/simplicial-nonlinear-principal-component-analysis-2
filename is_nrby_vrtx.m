function is_nrby_vrtx = is_nrby_vrtx(vrtx_crdnts, cntr_crdnts, dstnc_tol)

if isempty(vrtx_crdnts)

    is_nrby_vrtx = false(0);
    return;

end

crdnts_dmnsn = size(vrtx_crdnts, 1);
extents      = zeros(crdnts_dmnsn, 2);

%search sphere is centered at cntr_crdnts with radius dstnc_tol
%build extents for tight box that encloses search sphere
% for k=1:crdnts_dmnsn
% 
%     extents(k, :) = ...
%         [cntr_crdnts(k,1) - dstnc_tol, cntr_crdnts(k,1) + dstnc_tol];
% 
% end

extents(:, 2) = cntr_crdnts(:);
extents(:, 1) = cntr_crdnts(:);
extents(:, 1) = extents(:, 1) - dstnc_tol;
extents(:, 2) = extents(:, 2) + dstnc_tol;

%find points inside box, then determine which points in the box are in the
%sphere
is_nrby_vrtx  = in_box(vrtx_crdnts, extents);
% is_nrby_vrtx_inds = find(is_nrby_vrtx);

% dstncs = zeros(size(is_nrby_vrtx_inds));

dsplcmnts  = bsxfun(@minus, vrtx_crdnts(:, is_nrby_vrtx), cntr_crdnts(:));
dstncsSqrd = sum(dsplcmnts.*dsplcmnts, 1);

% for k=1:numel(dstncs)
% 
%     dstncs(k) = ...
%         norm(vrtx_crdnts(:, is_nrby_vrtx_inds(k)) - cntr_crdnts(:));
%     
% end
is_nrby_vrtx_inds = find(is_nrby_vrtx);
is_nrby_vrtx(is_nrby_vrtx_inds(dstncsSqrd > dstnc_tol*dstnc_tol)) = false;
