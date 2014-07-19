function [cand_tri_ovrlap, exstng_tri_ovrlap] ...
    = ...
    cand_tri_overlap_0_shrd_vrtcs(...
    cand_tri_vrtx_crdnts, ...
    vrtx_crdnts, ...
    exstng_tri_vrtx_inds)

num_exstng_tris = size(exstng_tri_vrtx_inds, 1);

cand_tri_ovrlap   = false(1, num_exstng_tris);
exstng_tri_ovrlap = false(1, num_exstng_tris);

cand_tri_vrtx_trnsltd_crdnts = [...
    cand_tri_vrtx_crdnts(:, 2) - cand_tri_vrtx_crdnts(:, 1), ...
    cand_tri_vrtx_crdnts(:, 3) - cand_tri_vrtx_crdnts(:, 1)];

cand_tri_plnr_vrtx_crdnts = zeros(2, 3);
[...
    tri_basis, ...
    cand_tri_plnr_vrtx_crdnts(:, 2:3)] ...
    = qr(cand_tri_vrtx_trnsltd_crdnts, 0);

for tri_i=1:num_exstng_tris

    %translate and project existing triangle into subspace containing 
    %candidate triangle
    exstng_tri_vrtx_trnsltd_crdnts = ...
        vrtx_crdnts(:, exstng_tri_vrtx_inds(tri_i, :));

    exstng_tri_vrtx_trnsltd_crdnts = bsxfun(...
        @minus, ...
        exstng_tri_vrtx_trnsltd_crdnts, ...
        cand_tri_vrtx_crdnts(:, 1));
    
            
    exstng_tri_plnr_vrtx_crdnts = ...
        tri_basis'*exstng_tri_vrtx_trnsltd_crdnts;
    
    cand_tri_ovrlap(tri_i) = plnr_tri_tri_cnflct(...
        cand_tri_plnr_vrtx_crdnts(1,:), ...
        cand_tri_plnr_vrtx_crdnts(2,:), ...
        exstng_tri_plnr_vrtx_crdnts(1,:), ...
        exstng_tri_plnr_vrtx_crdnts(2,:), ...
        [], []);

    %\/ PLOT DEBUG \/
    %     figure(2); subplot(1,2,1); plot_0_ce_shrd_hndl = plot(...
    %         cand_tri_plnr_vrtx_crdnts(1, [1:3, 1]), ...
    %         cand_tri_plnr_vrtx_crdnts(2, [1:3, 1]), ...
    %         'b-', ...
    %         exstng_tri_plnr_vrtx_crdnts(1, [1:3, 1]), ...
    %         exstng_tri_plnr_vrtx_crdnts(2, [1:3, 1]), ...
    %         'g-');
    %     axis equal;
    %     delete(plot_0_ce_shrd_hndl);
    %/\ PLOT DEBUG /\

    
end

exstng_tri_plnr_vrtx_crdnts  = zeros(2, 3);

for tri_i=1:num_exstng_tris

    %translate and project candidate triangle into subspace containing
    %existing triangle
    exstng_tri_vrtx_trnsltd_crdnts  = ...
        vrtx_crdnts(:, exstng_tri_vrtx_inds(tri_i, 2:3));
    
    exstng_tri_vrtx_trnsltd_crdnts = bsxfun(...
        @minus, ...
        exstng_tri_vrtx_trnsltd_crdnts, ...
        vrtx_crdnts(:, exstng_tri_vrtx_inds(tri_i, 1)));
    
    cand_tri_vrtx_trnsltd_crdnts = cand_tri_vrtx_crdnts;
    
    cand_tri_vrtx_trnsltd_crdnts = bsxfun(...
        @minus, ...
        cand_tri_vrtx_trnsltd_crdnts, ...
        vrtx_crdnts(:, exstng_tri_vrtx_inds(tri_i, 1)));
        
 
    [tri_basis, exstng_tri_plnr_vrtx_crdnts(:, 2:3)] ...
        = qr(exstng_tri_vrtx_trnsltd_crdnts, 0);
               
    cand_tri_plnr_vrtx_crdnts ...
        = tri_basis'*cand_tri_vrtx_trnsltd_crdnts;
        
    exstng_tri_ovrlap(tri_i) = plnr_tri_tri_cnflct(...
        cand_tri_plnr_vrtx_crdnts(1,:), ...
        cand_tri_plnr_vrtx_crdnts(2,:), ...
        exstng_tri_plnr_vrtx_crdnts(1,:), ...
        exstng_tri_plnr_vrtx_crdnts(2,:), ...
        [], []);
    
        %\/ PLOT DEBUG \/
        %     figure(2); subplot(1,2,1); plot_0_ce_shrd_hndl = plot(...
        %         cand_tri_plnr_vrtx_crdnts(1, [1:3, 1]), ...
        %         cand_tri_plnr_vrtx_crdnts(2, [1:3, 1]), ...
        %         'b-', ...
        %         exstng_tri_plnr_vrtx_crdnts(1, [1:3, 1]), ...
        %         exstng_tri_plnr_vrtx_crdnts(2, [1:3, 1]), ...
        %         'g-');
        %     axis equal;
        %     delete(plot_0_ce_shrd_hndl);
    %/\ PLOT DEBUG /\

    
end
