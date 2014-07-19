function [vrtxCrdnts, triVrtxInds, edgVrtxInds, ...
    P_dmnt_egnvctrs, P_dmnt_egnvals]...
    = SNPCA_interleaved_main(srfcCrdnts, SNPCA_params)

vrtxCrdnts  = zeros(size(srfcCrdnts, 1), SNPCA_params.INTL_NUM_VRTXS);
edgVrtxInds = zeros(SNPCA_params.INTL_NUM_EDGS, 3);
triVrtxInds = zeros(SNPCA_params.INTL_NUM_TRIS, 3);

numVrtxs = 0;
numEdgs  = 0;
numTris  = 0;

%emprcl_drctn_crrltn_mtrx{vrtx_ind} holds the empirical local direction
%correlation matrix at the vertex associated with vrtx_ind
%emprcl_drctn_crrltn_mtrx = cell(1, SNPCA_params.MAX_NUM_VRTCS);

num_srfc_data_pnts = size(srfcCrdnts, 2);

%empirical covariance matrix P = P_factor*P_factor'
%P_factor        = cell(1, SNPCA_params.MAX_NUM_VRTCS);
%use map with two keys (one key string from catenated ieee to hex function)
P_dmnt_egnvctrs = cell(SNPCA_params.INTL_NUM_EDGS, 2);
P_dmnt_egnvals  = cell(SNPCA_params.INTL_NUM_EDGS, 2);

%tell eigs that empirical covariance matrices are real symmetric
eigs_opts       = struct('issym', 1, 'isreal', 1);


%store minimizer computed by constrained optimization problem so we only
%compute it once
%the coordinates of the candidate vertex generated from vertices with
%indices v1 and v2 are stored at intl_mnmzr_crdnts{v1,v2} and
%intl_mnmzr_crdnts{v2,v1}
%intl_mnmzr_crdnts = cell(SNPCA_params.MAX_NUM_VRTCS); %space inefficient!
intl_mnmzr_crdnts = {};

edg_fifo = zeros(1, SNPCA_params.INTL_EDG_FIFO_LNGTH);

srfc_pnt_is_vbl_intl = true(size(srfcCrdnts, 2), 1);

%surf_pt_blngs_to_tri:
%does the surface point belong to a triangle (boolean)
%if so, what triangle (triangles?) does it belong to?
%when placing a triangle vertex candidate, we may not want to consider
%surface points that near the surface of an existing triangle
% surf_pt_blngs_to_tri(num_surf_pts) = ...
%     struct('belongs', false(num_surf_pts,1), 'tri_inds', []);

plot_hndls = new_plot_hndls();
if SNPCA_params.plot_frqncy > 0
    %plot surface points
    figure(1);
    [az, el] = view();
    clf
    hold off
    
    
    %surface points
    srfc_crdnts_3D = SNPCA_params.rtn_mtrx(:,1:3).'*srfcCrdnts;
    plot_hndls.surf_pts_hndl = plot_surf_pts(...
        srfc_crdnts_3D(1, :), ...
        srfc_crdnts_3D(2, :), ...
        srfc_crdnts_3D(3, :));
    %set(plot_hndls.surf_pts_hndl, 'Visible', 'off')
    clear srfc_crdnts_3D ;
    
    %warning('Surface data points are invisible')
    %set(plot_hndls.surf_pts_hndl, 'Visible', 'off');
    
    %disp('Clearing axis ticks')
    set(gca, ...
        'XTick', [], ...
        'YTick', [], ...
        'ZTick', [])
    
    %use open gl renderer for speed
    set(gcf, 'renderer', 'opengl')
    
    
    view(az, el);
    axis equal
    axis vis3d;
    
    hold on
end


%\/ Generate initial triangulation \/
%gen_init_tris2 generates an edge and two vertices so that the edge and the
%two vertices form the two initial triangles
[intl_edg_srfc_inds, intl_vrtx_inds] = gen_init_tris2(...
    srfcCrdnts, ...
    SNPCA_params.intl_pt_ind, ...
    SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.cnstrnt_rad_fac*SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.chrctrstc_lngth, SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.emprcl_drctn_crrltn_eval_bias, ...
    SNPCA_params.emprcl_drctn_egn_sprs_algrthm, ...
    SNPCA_params.nnz_egnvals);

numTris = numel(intl_vrtx_inds);

if numTris == 0
   
    error('Error generating initial triangulation.')
    
else
    
    vrtx_ind_to_srfc_pt_ind = ...
        [intl_edg_srfc_inds(:).' intl_vrtx_inds(:).'];
    
    srfc_pt_ind_to_vrtx_ind = ...
        containers.Map(num2cell(vrtx_ind_to_srfc_pt_ind), ...
                       num2cell(1:numel(vrtx_ind_to_srfc_pt_ind)), ...
                       'uniformValues', true);
    %works on newer versions of Matlab    
    %     srfc_pt_ind_to_vrtx_ind = ...
    %         containers.Map(vrtx_ind_to_srfc_pt_ind, ...
    %         1:numel(vrtx_ind_to_srfc_pt_ind));
    %works on newer versions of Matlab

    numVrtxs = 2 + numel(intl_vrtx_inds);
    
    vrtxCrdnts(:, 1:numVrtxs) = srfcCrdnts(:, vrtx_ind_to_srfc_pt_ind);

    switch numTris
        case 2            
            
            numEdgs = 5;
            edgVrtxInds(1:numEdgs,:) = [...
                1 2 3; ...
                2 3 1; ...
                1 3 2; ...
                2 4 1; ...
                1 4 2];

            triVrtxInds(1:2,:) = [...
                1 2 3; ...
                1 2 4];
            
            %push front edges of initial triangulation onto the edge stack
            edg_fifo(1:4) = [2 3 4 5];
            edg_fifo_ind  = 4;

        case 1
        
            numEdgs = 3;
            edgVrtxInds(1:numEdgs, :) = [...
                1 2 2; ...
                2 3 1; ...
                1 3 2];
            
            triVrtxInds(1, :) = [1 2 3];
        
        %push front edges of initial triangulation onto the edge stack
        edg_fifo(1:3) = [1 2 3];
        edg_fifo_ind  = 3;
        
            
    end
    
end

%/\ Generate initial triangulation /\


data_dmnsn = size(srfcCrdnts, 1);

%cndt_vrtx_info = new_cndt_vrtx_info(data_dmnsn);



if SNPCA_params.plot_frqncy > 0
    
    %plot triangulation
    
    vrtx_crdnts_3D = ...
        SNPCA_params.rtn_mtrx(:,1:3).'*vrtxCrdnts(:,1:numVrtxs);
    plot_hndls.tris_hndl = plot_tris(...
        triVrtxInds(1:numTris, :), ...
        vrtx_crdnts_3D(1,:), ...
        vrtx_crdnts_3D(2,:), ...
        vrtx_crdnts_3D(3,:), []);
    
end


%sav_cnt = 0;

%dbg_num_tris      = size(tri_vrtx_inds,1);
%dbg_prvs_num_tris = dbg_num_tris;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BEGIN load old program state and replot
%   %load; subplot(2,1,1); hold on; axis square; axis equal;
% load('saved_data/nt_2.mat'); %all_frnt_cycls_are_dead = false; seed is
% %load
%  plot_progress = true;
% %  load('tris_data_1700.mat'); plot_progress = true; clf; hold on; axis equal
% % seam_edg_inds = [];
% % front_is_subset_of_seams = false;
% surf_pts_hndl           = plot3([], [], [], '');
% actv_edge_hndl          = plot3([], [], [], '');
% actv_edge_intr_hndl     = plot3([], [], [], '');
% cand_vert_and_surf_hndl = plot3([], [], [], '');
% front_hndls             = plot3([], [], [], '');
% tris_hndl = plot_tris(tri_vrtx_inds, vrtx_crds_x, vrtx_crds_y, vrtx_crds_z);
%END load old program state and replot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% profile off
% profile on

%\/ skip by calling load \/
if true
    [...
        vrtxCrdnts, triVrtxInds, edgVrtxInds, ...
        numVrtxs, numTris, numEdgs, ...
        srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind,...
        intl_mnmzr_crdnts, ...
        P_dmnt_egnvctrs, P_dmnt_egnvals...
        ]...
        = ...
        advancing_front_main_loop(...
        srfcCrdnts, ...
        edg_fifo, edg_fifo_ind, ...
        vrtxCrdnts, triVrtxInds, edgVrtxInds, ...
        numVrtxs, numTris, numEdgs, ...
        srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind,...
        intl_mnmzr_crdnts, ...
        P_dmnt_egnvctrs, P_dmnt_egnvals, eigs_opts, ...
        SNPCA_params, plot_hndls);
    
    %mark all surface data points that coincide with vertices as inviable
    %starting points
    srfc_pnt_is_vbl_intl(vrtx_ind_to_srfc_pt_ind) = false;
    
    if SNPCA_params.save_data
        disp('Saving initial advancing front run.');
        save(...
            [SNPCA_params.path_saved_data ...
            'SNPCA_interleaved_run_' datestr(now, 30)]);
    end
    
else
    
    %load('interleaved_data/SNPCA_interleaved_run_20131011T111550.mat');
    load('interleaved_data/SNPCA_interleaved_run_20131024T134914.mat');    
    %SNPCA_params.non_adj_tri_ovrlap_dist_tol = .1*SNPCA_params.chrctrstc_lngth;
    %SNPCA_params.non_adj_tri_dist_tol = .025*SNPCA_params.chrctrstc_lngth;
    %SNPCA_params.new_tri_max_edg_lngth = 2*SNPCA_params.chrctrstc_lngth;

    vrtx_crdnts_3D = SNPCA_params.rtn_mtrx(:,1:3)'*vrtxCrdnts(:, 1:numVrtxs);
    plot_hndls.tris_hndl = plot_tris(triVrtxInds(1:numTris,:), ...
        vrtx_crdnts_3D(1,:), ...
        vrtx_crdnts_3D(2,:), ...
        vrtx_crdnts_3D(3,:), ...
        plot_hndls.tris_hndl);
    %set plot background color to white
    set(gcf, 'Color', [1 1 1]);

end

if true
    [triVrtxInds, edgVrtxInds, numTris, numEdgs] = sew_seams_decomp4(...
        srfcCrdnts, ...
        vrtx_ind_to_srfc_pt_ind, ...
        vrtxCrdnts, ...
        triVrtxInds, ...
        edgVrtxInds, ...
        numVrtxs, ...
        numTris, ...
        numEdgs, ...
        [], ...
        0, ...
        P_dmnt_egnvctrs, ...
        P_dmnt_egnvals, ...
        SNPCA_params.emprcl_drctn_crrltn_eval_bias, ...
        SNPCA_params.nnz_egnvals, ...
        SNPCA_params.emprcl_drctn_egn_sprs_algrthm, ...
        eigs_opts, ...
        SNPCA_params.new_tri_max_edg_lngth, ...
        SNPCA_params.non_adj_tri_dist_tol, ...
        SNPCA_params.non_adj_tri_ovrlap_dist_tol, ...
        SNPCA_params.rtn_mtrx, ...
        plot_hndls, ...
        SNPCA_params.plot_frqncy);
        
    if SNPCA_params.save_data
        
        disp('Saving initial seam sewing run.');
        disp('Saving initial advancing front run.');
        save(...
            [SNPCA_params.path_saved_data ...
            'SNPCA_interleaved_run_' datestr(now, 30)]);
    end
    
    if SNPCA_params.plot_frqncy > 0
        vrtx_crdnts_3D = ...
            SNPCA_params.rtn_mtrx(:,1:3)'*vrtxCrdnts(:, 1:numVrtxs);
        plot_hndls.tris_hndl = plot_tris(...
            triVrtxInds(1:numTris, :), ...
            vrtx_crdnts_3D(1,:), ...
            vrtx_crdnts_3D(2,:), ...
            vrtx_crdnts_3D(3,:), ...
            plot_hndls.tris_hndl);
    end
    
else
    
    %load('interleaved_data/SNPCA_interleaved_run_20130807T173239.mat');
    load('interleaved_data/SNPCA_interleaved_run_20131024T142408.mat');
    
    
    vrtx_crdnts_3D = SNPCA_params.rtn_mtrx(:,1:3)'*vrtxCrdnts(:, 1:numVrtxs);
    plot_hndls.tris_hndl = plot_tris(triVrtxInds(1:numTris, :), ...
        vrtx_crdnts_3D(1,1:numVrtxs), ...
        vrtx_crdnts_3D(2,1:numVrtxs), ...
        vrtx_crdnts_3D(3,1:numVrtxs), ...
        plot_hndls.tris_hndl);
    set(gcf, 'Color', [1 1 1]);
    
end

fnd_vbl_intl_srfc_pnt = true;
plot_hndls_to_delete = {...
    'actv_edge_hndl', 'actv_edge_intr_hndl', 'front_hndls', ...
    'cand_vert_and_surf_hndl'};

restart_count = 0;

while fnd_vbl_intl_srfc_pnt ...
        && restart_count <= (SNPCA_params.max_num_restarts-1);

    
    restart_count = restart_count + 1;
    if SNPCA_params.save_data
        
        disp(['Saving. Number of restarts: ' num2str(restart_count)]);
        save(...
            [SNPCA_params.path_saved_data ...
            'SNPCA_interleaved_run_' datestr(now, 30)]);

    end
    

    [intl_tri_is_vbl, ...
    intl_edg_srfc_inds, ...
    intl_vrtx_inds,...
    srfc_pnt_is_vbl_intl] ...
    ...
    = intl_srfc_data_pnt(...
    ...
    srfcCrdnts, ...
    vrtxCrdnts, ...
    triVrtxInds, ...
    edgVrtxInds, ...
    srfc_pnt_is_vbl_intl, ...
    SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.cnstrnt_rad_fac*SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.chrctrstc_lngth, ...
    SNPCA_params.emprcl_drctn_crrltn_eval_bias, ...
    SNPCA_params.nnz_egnvals, ...
    SNPCA_params.emprcl_drctn_egn_sprs_algrthm, ...
    SNPCA_params.non_adj_tri_dist_tol, ...
    SNPCA_params.nrby_vrtx_dstnc_tol, ...
    SNPCA_params.non_adj_tri_ovrlap_dist_tol, ...
    num_srfc_data_pnts);

    fnd_vbl_intl_srfc_pnt = any(intl_tri_is_vbl);


    %restart if there's a data point that's far from the triangulation and
    %viable
    if fnd_vbl_intl_srfc_pnt
               
        % \/ add viable initial triangles to the existing triangulation \/        
        [...
            new_tri_vrtx_inds, new_edg_vrtx_inds, newVrtxInds, ...
            new_vrtx_crdnts, ...
            srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind] ...
            = ...
            intl_cndt_tris_update(...
            intl_tri_is_vbl, intl_edg_srfc_inds, intl_vrtx_inds, ...
            numVrtxs, srfcCrdnts, ...
            srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind);

        numNewVrtxs = numel(newVrtxInds);
        numNewEdgs  = size(new_edg_vrtx_inds, 1);
        numNewTris  = size(new_tri_vrtx_inds, 1);
        
        triVrtxInds((numTris+1):(numTris+numNewTris), :) ...
            = new_tri_vrtx_inds;
        edgVrtxInds((numEdgs+1):(numEdgs+numNewEdgs), :) ...
            = new_edg_vrtx_inds;
        vrtxCrdnts(:, (numVrtxs+1):(numVrtxs+numNewVrtxs)) ...
            = new_vrtx_crdnts;        
        
        numVrtxs = numVrtxs + numNewVrtxs;
        numEdgs  = numEdgs + numNewEdgs;
        numTris  = numTris + numNewTris;
        
        edg_fifo(1:end) = 0;
        edg_fifo(1:numNewEdgs) = (numEdgs-numNewEdgs+1):numEdgs;
        edg_fifo_ind  = numNewEdgs;
        % /\ add viable initial triangles to the existing triangulation /\

        
        [...
            vrtxCrdnts, triVrtxInds, edgVrtxInds, ...
            numVrtxs, numTris, numEdgs, ...
            srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind,...
            intl_mnmzr_crdnts, ...
            P_dmnt_egnvctrs, P_dmnt_egnvals...
            ]...
            = ...
            advancing_front_main_loop(...
            srfcCrdnts, ...
            edg_fifo, edg_fifo_ind, ...
            vrtxCrdnts, triVrtxInds, edgVrtxInds, ...
            numVrtxs, numTris, numEdgs, ...
            srfc_pt_ind_to_vrtx_ind, vrtx_ind_to_srfc_pt_ind,...
            intl_mnmzr_crdnts, ...
            P_dmnt_egnvctrs, P_dmnt_egnvals, eigs_opts, ...
            SNPCA_params, plot_hndls);

        
        %delete plot detritus
        if exist('plot_hndls', 'var') && isstruct(plot_hndls)
            for k=1:numel(plot_hndls_to_delete)
                
                if isfield(plot_hndls, plot_hndls_to_delete{k}) ...
                        && ~isempty(plot_hndls.(plot_hndls_to_delete{k})) ...
                        && ishandle(plot_hndls.(plot_hndls_to_delete{k}))
                    %plot_hndls.(plot_hndls_to_delete{k}) "dynamic field",
                    %access by string name of field
                    delete(plot_hndls.(plot_hndls_to_delete{k}));
                end
                
            end
            
        end

        [triVrtxInds, edgVrtxInds, numTris, numEdgs] = sew_seams_decomp4(...
            srfcCrdnts, ...
            vrtx_ind_to_srfc_pt_ind, ...
            vrtxCrdnts, ...
            triVrtxInds, ...
            edgVrtxInds, ...
            numVrtxs, ...
            numTris, ...
            numEdgs, ...
            [], ...
            0, ...
            P_dmnt_egnvctrs, ...
            P_dmnt_egnvals, ...
            SNPCA_params.emprcl_drctn_crrltn_eval_bias, ...
            SNPCA_params.nnz_egnvals, ...
            SNPCA_params.emprcl_drctn_egn_sprs_algrthm, ...
            eigs_opts, ...
            SNPCA_params.new_tri_max_edg_lngth, ...
            SNPCA_params.non_adj_tri_dist_tol, ...
            SNPCA_params.non_adj_tri_ovrlap_dist_tol, ...
            SNPCA_params.rtn_mtrx, ...
            plot_hndls, ...
            SNPCA_params.plot_frqncy);
        

        %delete plot detritus
        if exist('plot_hndls', 'var') && isstruct(plot_hndls)
            for k=1:numel(plot_hndls_to_delete)
                
                if isfield(plot_hndls, plot_hndls_to_delete{k}) ...
                        && ~isempty(plot_hndls.(plot_hndls_to_delete{k})) ...
                        && ishandle(plot_hndls.(plot_hndls_to_delete{k}))
                    %plot_hndls.(plot_hndls_to_delete{k}) "dynamic field",
                    %access by string name of field
                    delete(plot_hndls.(plot_hndls_to_delete{k}));
                end
                
                
            end
            
        end
        
        if SNPCA_params.plot_frqncy > 0
            vrtx_crdnts_3D = SNPCA_params.rtn_mtrx(:,1:3)'*vrtxCrdnts(:, 1:numVrtxs);
            plot_hndls.tris_hndl = plot_tris(...
                triVrtxInds(1:numTris, :), ...
                vrtx_crdnts_3D(1,:), ...
                vrtx_crdnts_3D(2,:), ...
                vrtx_crdnts_3D(3,:), ...
                plot_hndls.tris_hndl );
        end
        %delete_plot_lines(gca);
        
    end
    
end


%draw the tessellation and front one last time

if SNPCA_params.plot_frqncy > 0
    
    vrtx_crdnts_3D = ...
        SNPCA_params.rtn_mtrx(:,1:3)'*vrtxCrdnts(:, 1:numVrtxs);
    plot_hndls.tris_hndl = plot_tris(triVrtxInds(1:numTris, :), ...
        vrtx_crdnts_3D(1,:), ...
        vrtx_crdnts_3D(2,:), ...
        vrtx_crdnts_3D(3,:), ...
        plot_hndls.tris_hndl);
    %set plot background color to white
    set(gcf, 'Color', [1 1 1])
    
end



% vrtx_crdnts_3D = NLPCA_params.rtn_mtrx(:,1:3).'*vrtx_crdnts;
%
% update_tri_surf(plot_hndls.tris_hndl, tri_vrtx_inds,...
%     vrtx_crdnts_3D(1,:), ...
%     vrtx_crdnts_3D(2,:), ...
%     vrtx_crdnts_3D(3,:));
%
% delete(plot_hndls.actv_edge_hndl(ishandle(plot_hndls.actv_edge_hndl)));
% delete(plot_hndls.actv_edge_intr_hndl(ishandle(...
%     plot_hndls.actv_edge_intr_hndl)));
% delete(plot_hndls.front_hndls(plot_hndls.ishandle(...
%     front_hndls)));
% delete(plot_hndls.cand_vert_and_surf_hndl(ishandle(...
%     plot_hndls.cand_vert_and_surf_hndl)));

% front_hndls = ...
%     plot_front(frnt_cycl_edg_inds, frnt_cycl_ind, ...
%     EdgVrtxInds, tri_vrtx_inds, front_info, ...
%     vrtx_crdnts_3D(1,:), ...
%     vrtx_crdnts_3D(2,:), ...
%     vrtx_crdnts_3D(3,:));


%save(['advancing_front_data' filesep() 'SNPCA_run_' datestr(now, 30)])

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%advancing_front_main_loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dstncs = pnt_vrtx_dstncs(...
    pnt_crdnts, ...
    vrtx_inds, ...
    vrtx_crdnts)

num_vrtxs = numel(vrtx_inds);
dstncs    = zeros(num_vrtxs,1);

if isempty(dstncs)
    return
end

dim_pnt_crdnts  = size(pnt_crdnts, 1);
dim_vrtx_crdnts = size(vrtx_crdnts, 1);

assert (dim_pnt_crdnts == dim_vrtx_crdnts);

dstncs = ( vrtx_crdnts(1, vrtx_inds) - pnt_crdnts(1)).^2;
for k=2:dim_vrtx_crdnts
    
    dstncs = dstncs + ( vrtx_crdnts(k, vrtx_inds) - pnt_crdnts(k) ).^2;
    
end

dstncs = sqrt(dstncs);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%a front edge is near pnt_crdnts if a vertex of the front edge is within
%dstnc_tol of pnt_crdnts
function nrby_frnt_edg_inds = nrby_frnt_edgs(...
    pnt_crdnts, edg_vrtx_inds, tri_vrtx_inds, vrtx_crdnts, dstnc_tol)


%build list of vertices that belong to a front edge
frnt_edg_inds = all_frnt_edg_inds(edg_vrtx_inds, tri_vrtx_inds);
frnt_edg_vrtx_inds = edg_vrtx_inds(frnt_edg_inds, 1:2);

unq_frnt_edg_vrtx_inds = unique(frnt_edg_vrtx_inds(:));

%compute distance from the point with coordinates pnt_crdnts to front edge
%vertices
dstncs = pnt_vrtx_dstncs(pnt_crdnts, unq_frnt_edg_vrtx_inds, vrtx_crdnts);

dstnc_lt_tlrnc = dstncs < dstnc_tol;

nrby_frnt_vrtx_inds = unq_frnt_edg_vrtx_inds(dstnc_lt_tlrnc);

num_nrby_frnt_vrtxs = numel(nrby_frnt_vrtx_inds);

if num_nrby_frnt_vrtxs == 0
    
    nrby_frnt_edg_inds = [];
    
end

is_nrby_frnt_edg = ...
    nrby_frnt_vrtx_inds(1) == frnt_edg_vrtx_inds(:, 1) ...
    | nrby_frnt_vrtx_inds(1) == frnt_edg_vrtx_inds(:, 2);


for k=2:num_nrby_frnt_vrtxs
    
    is_nrby_frnt_edg = ...
        is_nrby_frnt_edg ...
        | ...
        nrby_frnt_vrtx_inds(k) == frnt_edg_vrtx_inds(:, 1) ...
        | nrby_frnt_vrtx_inds(k) == frnt_edg_vrtx_inds(:, 2);
    
end

nrby_frnt_edg_inds = frnt_edg_inds(is_nrby_frnt_edg);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


