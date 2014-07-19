function QMtrxVctrPrdct = QMtrxVctrPrdctDecomp(...
    x, Q_small_egnvctrs, Q_small_egnvals, Q_large_egnval)

x_crndt_vctr = Q_small_egnvctrs'*x(:);
Q_small_egnvctrs_x_crndt_vctr_prdct = Q_small_egnvals(:).*x_crndt_vctr;
QMtrxVctrPrdct = Q_small_egnvctrs*Q_small_egnvctrs_x_crndt_vctr_prdct;

[numAllEgnvctrs, numSmallEgnvctrs] = size(Q_small_egnvctrs);
if numAllEgnvctrs > numSmallEgnvctrs
    %numAllEgnvctrs > numSmallEgnvctrs test avoids arithmetic with
    %Q_large_egnval, i.e. when Q_large_egnval = Inf

    x_perp = x(:) - Q_small_egnvctrs*(Q_small_egnvctrs*x_crndt_vctr);
    
    QMtrxVctrPrdct = QMtrxVctrPrdct + Q_large_egnval*x_perp;

end
