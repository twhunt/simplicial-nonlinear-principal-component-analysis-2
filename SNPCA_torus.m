function [] = SNPCA_torus()

state_space_dmnsn = 3;  

SNPCA_params = new_SNPCA_params();

SNPCA_params.chrctrstc_lngth               = .75; 
SNPCA_params.cnstrnt_rad_fac               = .5*sqrt(3);
SNPCA_params.new_tri_max_edg_lngth         = 2*SNPCA_params.chrctrstc_lngth;
SNPCA_params.srch_rad_fac1                 = 1.0;
SNPCA_params.srch_rad_fac2                 = 1.0;
SNPCA_params.intl_pt_ind                   = 1;
SNPCA_params.prfrd_cnstrnt_rds             = .5*sqrt(3)*SNPCA_params.chrctrstc_lngth;
SNPCA_params.prfrd_cnstrnt_rds_wght        = .5;
SNPCA_params.non_adj_tri_dist_tol          = .2*SNPCA_params.chrctrstc_lngth;
SNPCA_params.non_adj_tri_ovrlap_dist_tol   = .25*SNPCA_params.chrctrstc_lngth;
SNPCA_params.nrby_vrtx_dstnc_tol           = 2.5*SNPCA_params.chrctrstc_lngth;
SNPCA_params.cand_vert_max_nudge_dist      = .25*SNPCA_params.chrctrstc_lngth;
SNPCA_params.adj_vert_max_nudge_dist       = .25*SNPCA_params.chrctrstc_lngth;    
SNPCA_params.emprcl_drctn_crrltn_eval_bias = 1e-6;
SNPCA_params.plot_frqncy                   = 1;
SNPCA_params.rtn_mtrx                      = eye(3);
%[SNPCA_params.rtn_mtrx, tmp_R]             = qr(rand(state_space_dmnsn));
SNPCA_params.nnz_egnvals                   = 3;
SNPCA_params.emprcl_drctn_egn_sprs_algrthm = false;
SNPCA_params.max_num_restarts              = 0;
SNPCA_params.save_data                     = true;


%Generate points on torus surface

%Set add_noise to true to add noise to the surface data coordinates
%Set add_noise_to_all_coordinates to true to add noise to each coordinate
%after rotation into the higher dimensional space. 
%Set it to false to add noise to the torus coordinates before rotaion into
%the higher dimensional space.
add_noise = false;
add_noise_to_all_coordinates = false;
if add_noise
    %noise_magnitude = 0;
    noise_magnitude = .05;
    %noise_magnitude = .1;
    %noise_magnitude = .0001;
else
    noise_magnitude = 0;
end

seed = 5;
%rnd_strm = RandStream('mt19937ar','seed', seed);
%RandStream.setGlobalStream(rnd_strm);
n   = 1e4; %number of data 
rh  = 4; %horizontal torus radius
rd  = 1; %vertical torus radius
trs_crdnts = gen_torus_data_pts2(rh, rd, n);

if add_noise && ~add_noise_to_all_coordinates
    noise = randn(3, n);
    noise = noise_magnitude*noise;
    
    trs_crdnts = trs_crdnts + noise;
end

srfc_crdnts = SNPCA_params.rtn_mtrx(:,1:3)*trs_crdnts;

if add_noise && add_noise_to_all_coordinates
    noise = randn(state_space_dmnsn, n);
    noise = noise_magnitude*noise;
    
    srfc_crdnts = srfc_crdnts + noise;
end


disp(['Torus large radius: ' num2str(rh)])
disp(['Torus small radius: ' num2str(rd)])
disp(['state space dimension: ' num2str(state_space_dmnsn)])


disp(['Number surface data points: ' num2str(n)])
disp(['Noise magnitude: ' num2str(noise_magnitude)])


%\/ temp \/
%[Utmp Stmp Vtmp] = svd(srfc_crdnts);
%/\ temp /\

[vrtx_crdnts, tri_vrtx_inds, edg_vrtx_inds, ...
    P_dmnt_egnvctrs, P_dmnt_egnvals] ... 
    = ...
    SNPCA_interleaved_main(srfc_crdnts, SNPCA_params);

