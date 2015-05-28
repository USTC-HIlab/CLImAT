function [corrected_rd, m_all_map, m_bin_map] = MAP_Correction(data_rd, data_map)
% 21/05/2013 by yzh
% Mappability correction for read counts

% subplot(2,1,1);
% plot(data_map,data_rd,'b.','MarkerSize',3);
corrected_rd = data_rd;
m_all_map = median(data_rd);
int_map = floor(data_map*100);
int_map_u = unique(int_map);
m_map = zeros(size(int_map_u));

for i = 1:length(int_map_u)
    tv = int_map == int_map_u(i);
    m_map(i) = median(data_rd(tv));
    corrected_rd(tv) = round(data_rd(tv)*m_all_map/(m_map(i)+eps));
end
if size(m_map, 1) > size(m_map, 2)
    int_map_u = int_map_u';
    m_map = m_map';
end
m_bin_map = [int_map_u; m_map];
% subplot(2,1,2);
% plot(data_map,corrected_rd,'b.','MarkerSize',3);