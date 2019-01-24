function Image_Large      = Image_Block_Expand(Image, NskipX, NskipY)
% Image_Block_Expand makes larger image of Image by replicating each pixel
% with a block of NskipX x NskipY pixels

[NY,NX]     = size(Image);

Image_Large    = zeros( NY*NskipY , NX*NskipX);
blockM         = ones(NskipX,NskipY);

for ix = 1:(NX-1)
    for iy = 1:(NY-1),
        ix_start = (ix-1) * NskipX + 1;
        ix_end   = (ix)   * NskipX ;
        iy_start = (iy-1) * NskipY + 1;
        iy_end   = (iy)   * NskipY ;
        
        Image_Large(iy_start:iy_end , ix_start:ix_end) = blockM * Image(iy,ix);
        
    end
end



end

