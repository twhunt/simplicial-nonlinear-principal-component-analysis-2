function [P_dmnt_egnvctrs, P_dmnt_egnvals] = pnts_to_egn_dcmp_dns(...
    cntr_crdnts, lcl_srfc_crdnts, nnz_egnvals)

num_pts_in_ngbrhd = size(lcl_srfc_crdnts, 2)

unitDirectionsMtrx = unitDirections(cntr_crdnts, lcl_srfc_crdnts);

if size(unitDirectionsMtrx, 2) >= nnz_egnvals
    
    [P_dmnt_egnvctrs, P_dmnt_egnvals, throwaway] ...
        = svd(unitDirectionsMtrx, 'econ');
    
    P_dmnt_egnvals = diag(P_dmnt_egnvals);

else
    %number of surface data points in search sphere less than number of
    %empirical directian covariance matrix nonzero eigenvalues
    %pad vector of eigenvalues with zeros
    [P_dmnt_egnvctrs, P_dmnt_egnvals, throwaway] ...
        = svd(unitDirectionsMtrx);
    
    if size(P_dmnt_egnvals, 2) ~= 1
        %if size(unitDirectionsMtrx, 2) is 1, then diag(P_dmnt_egnvals)
        %returns a matrix instead of a vector
        P_dmnt_egnvals = diag(P_dmnt_egnvals);
    end

    P_dmnt_egnvals(end+1:nnz_egnvals) = 0;
    
end


P_dmnt_egnvals ...
    = (1/num_pts_in_ngbrhd)*(P_dmnt_egnvals.*P_dmnt_egnvals);

P_dmnt_egnvals  = P_dmnt_egnvals(1:nnz_egnvals);
P_dmnt_egnvctrs = P_dmnt_egnvctrs(:, 1:nnz_egnvals);

return

%\/ Not very wise version that forms the covariance matrix then
%diagonalizes instead of just computing the svd of the data matrix
emprcl_drctn_cvrnc_mtrx = ...
    emprcl_drctn_cvrnc(cntr_crdnts, lcl_srfc_crdnts);


% [P_dmnt_egnvctrs, tmp_diag_evals] = ...
%     eig(emprcl_drctn_cvrnc_mtrx, 'nobalance');

%emprcl_drctn_cvrnc_mtrx is symmetric positive semidefinite, but eig may
%return complex eigenvectors when emprcl_drctn_cvrnc_mtrx is singular (it
%happens in practice)
%This can't happen with svd, and right singular vectors are eigenvectors in
%the positive semidefinite case, so use svd.

[P_dmnt_egnvctrs_unused, tmp_diag_evals, P_dmnt_egnvctrs] = ...
    svd(emprcl_drctn_cvrnc_mtrx);

%Put eigenvalues in a 1D vector (instead of the diagonal matrix returned by
%svd)
P_dmnt_egnvals = diag(tmp_diag_evals);
%[P_dmnt_egnvals sort_prmtn] = sort(P_dmnt_egnvals, 'descend');

%only keep the dominant eigenvalues and eigenvectors specified by
%nnz_egnvals
P_dmnt_egnvals  = P_dmnt_egnvals(1:nnz_egnvals);
P_dmnt_egnvctrs = P_dmnt_egnvctrs(:, 1:nnz_egnvals);
% %eig called with nobalance flag produces eigenvectors with non-unit length
% for k=1:size(P_dmnt_egnvctrs)
%     
%     P_dmnt_egnvctrs(:,k) = ...
%         (1/norm(P_dmnt_egnvctrs(:,k)))*P_dmnt_egnvctrs(:,k);
%     
% end

