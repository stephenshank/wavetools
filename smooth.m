function S=smooth(M,smoothness)

S=imfilter(M,fspecial('gaussian',[50 50],smoothness),'replicate');