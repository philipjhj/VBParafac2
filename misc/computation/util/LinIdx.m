function ndx = LinIdx(dims,sub_idx)
ndx = sub_idx(1);
ndx = ndx + (sub_idx(2) - 1).*dims(1);
k = cumprod(dims);
for i = 3:numel(sub_idx)
    v = sub_idx(i);
    
    ndx = ndx + (v-1)*k(i-1);
end
end