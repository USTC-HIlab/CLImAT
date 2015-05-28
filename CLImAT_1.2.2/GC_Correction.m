function [corrected_rd, m_all_gc, m_bin_gc] = GC_Correction(data_rd, data_gc)
% 21/05/2013 by yzh
% GC correction for read counts

% subplot(2,1,1);
% plot(data_gc,data_rd,'b.','MarkerSize',4);
corrected_rd = data_rd;
m_all_gc = median(data_rd);
int_gc = floor(data_gc*100);
int_gc_u = unique(int_gc);
m_gc = zeros(size(int_gc_u));

for i = 1:length(int_gc_u)
    tv = int_gc == int_gc_u(i);
    m_gc(i) = median(data_rd(tv));
    corrected_rd(tv) = round(data_rd(tv)*m_all_gc/(m_gc(i)+eps));
end
if size(m_gc, 1) > size(m_gc, 2)
    int_gc_u = int_gc_u';
    m_gc = m_gc';
end
m_bin_gc = [int_gc_u; m_gc];
% subplot(2,1,2);
% plot(data_gc,corrected_rd,'b.','MarkerSize',4);