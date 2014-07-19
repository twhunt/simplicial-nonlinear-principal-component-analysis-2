function [sphrCnstrntVal, sphrCnstrntGrad] = sphereConstraint(x, rds)

sphrCnstrntVal  = sqrt(sum(x.^2));
sphrCnstrntGrad = (1/(rds*sphrCnstrntVal))*x(:);
sphrCnstrntVal  = sphrCnstrntVal/rds - 1;
