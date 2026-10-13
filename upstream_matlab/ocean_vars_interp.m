
function [temp,salt,pres,Ctemp,Sa] = ocean_vars_interp_v2(salt0,temp0,pres0,target,lon,lat)
n = size(salt0);
len = n(1,2);
wid = length(target);
temp = nan(wid,len);
salt = nan(wid,len);
parfor i = 1:len
    s = salt0(:,i);
    t = temp0(:,i);
    p = pres0(:,i);
    lonn = lon(i);
    latt = lat(i);
    salt1 = nan(wid,1);
    temp1 = nan(wid,1);
    [s1,t1,p1] = ocean_constant_pres_vars(s,t,p);
    if sum(~isnan(s1))>3
        % salt(:,i) = interp1(p1,s1,target);
        [sa,~] = gsw_SA_from_SP(s1,p1,lonn,latt);
        ct = gsw_CT_from_t(sa,t1,p1);
        [salt1,~] = gsw_SA_CT_interp(sa,ct,p1,target);
        salt1 = gsw_SP_from_SA(salt1,target,lonn,latt);
        idx1 = find(min(abs(p1(1)-target))==abs(p1(1)-target));
        idx2 = find(min(abs(p1(end)-target))==abs(p1(end)-target));
        salt1(1:idx1) = nan;
        salt1(idx2:end) = nan;
    elseif sum(~isnan(s1))<=3&&sum(~isnan(s1))>1
        [sa,~] = gsw_SA_from_SP(s,p,lonn,latt);
        ct = gsw_CT_from_t(sa,t,p);
        [salt1,~] = gsw_SA_CT_interp(sa,ct,p,target);
        salt1 = gsw_SP_from_SA(salt1,target,lonn,latt);
        idx1 = find(min(abs(p1(1)-target))==abs(p1(1)-target));
        idx2 = find(min(abs(p1(end)-target))==abs(p1(end)-target));
        salt1(1:idx1) = nan;
        salt1(idx2:end) = nan;
    elseif sum(~isnan(s1))==1
        [vs,idxs] = find_nearest(s1,p1,target);
        % salt(idxs,i) = vs;
        salt1(idxs) = vs;
    end
    if sum(~isnan(t1))>3
        % temp(:,i) = interp1(p1,t1,target);
        if sum(~isnan(s1))~=0
            [sa,~] = gsw_SA_from_SP(s1,p1,lonn,latt);
            ct = gsw_CT_from_t(sa,t1,p1);
            [~,temp1] = gsw_SA_CT_interp(sa,ct,p1,target);
            temp1 = gsw_t_from_CT(salt1,temp1,target);
        else
            try
                temp1 = gsw_t_interp(t1,p1,target);
            catch
                temp1 = gsw_t_interp(t,p,target);
            end
        end
    elseif sum(~isnan(t1))<=3&&sum(~isnan(t1))>1
        if sum(~isnan(s1))~=0
            [sa,~] = gsw_SA_from_SP(s,p,lonn,latt);
            ct = gsw_CT_from_t(sa,t,p);
            [~,temp1] = gsw_SA_CT_interp(sa,ct,p,target);
            temp1 = gsw_t_from_CT(salt1,temp1,target);
            
        else
            try
                temp1 = gsw_t_interp(t1,p1,target);
            catch
                temp1 = gsw_t_interp(t,p,target);
            end
        end
        idx1 = find(min(abs(p1(1)-target))==abs(p1(1)-target));
        idx2 = find(min(abs(p1(end)-target))==abs(p1(end)-target));
        temp1(1:idx1) = nan;
        temp1(idx2:end) = nan;
    elseif sum(~isnan(t1))==1
        [vt,idxt] = find_nearest(t1,p1,target);
        % temp(idxt,i) = vt;
        temp1(idxt) = vt;
    end
    temp(:,i) = temp1;
    salt(:,i) = salt1;
end
pres = repmat(target, 1, len);
Sa = gsw_SA_from_SP(salt,pres,lon,lat);
Ctemp = gsw_CT_from_t(Sa,temp,pres);
end