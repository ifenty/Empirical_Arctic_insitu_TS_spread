function [s1,t1,p1] = ocean_constant_pres_vars(s,t,p)
idx_p = ~isnan(p);
s1 = s(idx_p);
t1 = t(idx_p);
p1 = p(idx_p);
p1(p1>=6000) = nan;
end