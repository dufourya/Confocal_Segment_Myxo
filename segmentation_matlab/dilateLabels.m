function L = dilateLabels(L,im,s,t)

pL = L;
d = Inf;

while d>100
    [D, idx] = bwdist(L);
    idx = L(idx);
    idx(D>1 | im<=t) = 0;
    L = idx;
    L(~imopen(L>0,ones(s,s,s))) = 0;    
    d = sum((L(:) - pL(:))~=0);
    pL = L;
end

[~, idx] = bwdist(L);
idx = L(idx);
L = imfill(L,'holes');
idx(~L) = 0;
L = idx;