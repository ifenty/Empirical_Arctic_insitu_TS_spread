function dens_qc = ocean_dens_quality_check(dens,pres)
n = size(dens);
wid = n(1,1);
len = n(1,2);
dens_qc = zeros(wid,len);
dens_qc(isnan(dens))=nan;
parfor i = 1:len
    d = dens(:,i);
    local_qc = zeros(wid, 1);
    if sum(~isnan(d))~=1
        p = pres(:,i);
        Gd = diff(d(~isnan(d))) ./ diff(p(~isnan(d)));
        pp = p(~isnan(d));
        pp = pp(2:end);
        ind = pp <= 30;
        local_qc(ind) = Gd(ind) > -0.03;
        ind = pp > 30 & pp <= 400;
        local_qc(ind) = Gd(ind) > -0.02;
        ind = pp > 400;
        local_qc(ind) = Gd(ind) > -0.01;
    else
        ind = find(~isnan(d));
        local_qc(ind) = local_qc(ind) + 1;
    end
    dens_qc(:,i) = local_qc;
end
end