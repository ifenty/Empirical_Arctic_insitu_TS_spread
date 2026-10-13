% function pres_qc = ocean_pres_quality_check(pres)
% n = size(pres);
% wid = n(1,1);
% len = n(1,2);
% pres_qc = zeros(wid,len);
% pres_qc(isnan(pres))=nan;
% parfor i = 1:len
% p = pres(:,i);
% local_qc = zeros(wid, 1);
% if sum(~isnan(p))~=1
% ind1 = diff(p(~isnan(p)))>0;
% ind2 = ~all(isnan(p));
% ind = ind1&ind2;
% local_qc(ind) = local_qc(ind)+1;
% else
% ind = find(~isnan(p));
% local_qc(ind) = local_qc(ind) + 1;
% end
% pres_qc(:,i) = local_qc;
% end

%
%
%
%
% function pres_qc = ocean_pres_quality_check(pres)
% n = size(pres);
% wid = n(1,1);
% len = n(1,2);
% pres_qc = zeros(wid,len);
% pres_qc(isnan(pres)) = nan;
% parfor i = 1:len
% p = pres(:,i);
% local_qc = zeros(wid, 1);
% pp = p(~isnan(p));
% max_value = pp(1);
% if sum(~isnan(p))~=1
%     for j = 2:length(pp)
%         if pp(j) < max_value  % 发现异常值
%             continue
%         else
%             max_value = pp(j);  % 更新最大值
%             local_qc(j) = local_qc(j)+1;
%         end
%     end
% else
%     ind = find(~isnan(p));
%     local_qc(ind) = local_qc(ind) + 1;
% end
% pres_qc(:,i) = local_qc;
% end
% end



function pres_qc = ocean_pres_quality_check(pres)
[wid, len] = size(pres);
pres_qc = zeros(wid, len);
pres_qc(isnan(pres)) = nan;  % 保持 NaN 值不变

parfor i = 1:len
    p = pres(:, i);  % 获取第 i 列
    local_qc = zeros(wid, 1);
    valid_idx = find(~isnan(p));  % 获取非 NaN 值的索引
    pp = p(valid_idx);  % 提取非 NaN 数据
    
    if isempty(pp)
        continue;  % 当前列全是 NaN，跳过
    end
    if pp(1)>=0
    local_qc(valid_idx(1)) = local_qc(valid_idx(1))+1;
    end
    max_value = pp(1);
    for j = 2:length(pp)
        if pp(j)>=0
%             local_qc(valid_idx(j-1)) = local_qc(valid_idx(j-1))+1;
            if pp(j) <= max_value
                if j==2
                    local_qc(valid_idx(j-1)) = local_qc(valid_idx(j-1))+1;
                else
                    continue
                end
            else
                max_value = pp(j);  % 更新最大值
                local_qc(valid_idx(j)) = local_qc(valid_idx(j))+1;  % 标记正常增长点
            end
        end
        
        % 如果只有一个有效数据点，标记为正常
        if length(pp) == 1
            local_qc(valid_idx) = 1;
        end
    end
    pres_qc(:, i) = local_qc;
end
end

