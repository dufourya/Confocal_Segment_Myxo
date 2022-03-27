function [me1, me2, me3] = hessianFilter3D(I, sigma)

me1 = zeros(size(I));
me2 = zeros(size(I));
me3 = zeros(size(I));

for k = 1:numel(sigma)
      
    F = imgaussfilt3(I,sigma(k));

    [Dx,Dy,Dz]=gradient(F);
    
    [Hxx, Hxy, Hxz] = gradient(Dx);
    [~,   Hyy, Hyz] = gradient(Dy);
    [~,     ~, Hzz] = gradient(Dz);
    
    [e1,e2,e3]=eig3volume(Hxx,Hxy,Hxz,Hyy,Hyz,Hzz);
    
    idx = abs(me1)<abs(e1)*sigma(k);
    me1(idx) = e1(idx)*sigma(k);
    
    idx = abs(me2)<abs(e2)*sigma(k);
    me2(idx) = e2(idx)*sigma(k);
    
    idx = abs(me3)<abs(e3)*sigma(k);
    me3(idx) = e3(idx)*sigma(k);
    
end