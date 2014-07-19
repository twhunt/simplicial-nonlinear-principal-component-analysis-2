%P = emprcl_drctn_cvrnc(cntr_crds, nghbrhd_pts_crds)
%cntr_crds is a 1D array
%nghbrhd_pts_crds is a matrix of coordinates where each column is a
%coordinate vector
function unitDirectionsMtrx = unitDirections(cntrCrdnts, nghbrhdPntsCrdnts)

if isempty(cntrCrdnts)
    error('NLPCA:Empty_empty_crds', 'Empty coordinate vector');
end

num_nghbrhd_pts = size(nghbrhdPntsCrdnts, 2);

if num_nghbrhd_pts == 0
    error(...
        'NLPCA:no_surf_pts_in_sphr', ...
        'No surface points in search sphere');
end

%get direction by subtracting center from each neighborhood point
unitDirectionsMtrx = bsxfun(@minus, nghbrhdPntsCrdnts, cntrCrdnts(:));

%compute inverse norm of each direction
columnInvrsNorms = 1./sqrt(sum(unitDirectionsMtrx.*unitDirectionsMtrx, 1));

%scale each direction to have unit length
unitDirectionsMtrx = bsxfun(@times, unitDirectionsMtrx, columnInvrsNorms);