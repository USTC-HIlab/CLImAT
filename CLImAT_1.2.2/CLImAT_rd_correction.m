function corrected_rd = CLImAT_rd_correction(data_rd, data_gc, data_map)
% 13/01/2014 by yzh
% GC correction for read counts

int_gc = floor(data_gc*100);
int_map = floor(data_map*100);
corrected_rd = data_rd;
m_all = median(data_rd);
for i = min(int_gc):max(int_gc)
    for j = min(int_map):max(int_map)
        tv = int_gc == i & int_map == j;
        if sum(tv) > 0
            m = median(data_rd(tv));
            corrected_rd(tv) = round(data_rd(tv)*m_all/(m+eps));
        end
    end
end

% corrected_rd = GC_Map_Correction(data_rd, int_gc, int_map);

end