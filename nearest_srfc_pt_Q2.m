function [minQDstnc, nrstInd] ...
    = nearest_srfc_pt_Q2(srfcCrdnts, dstncMtrcDecomp)


numSrfcPnts = size(srfcCrdnts, 2);
for k=numSrfcPnts:-1:1
   
    Q_dstncs(k) = dstncMtrcDecomp(srfcCrdnts(:, k));        
    
end

[minQDstnc, nrstInd] = min(Q_dstncs);
