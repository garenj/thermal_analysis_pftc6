function [im_scaled] = rescale_image(im, qlo, qhi)
    %imss = im(~isnan(im));
    imss = im(im>0); % chop everything below -10... this may need to change for other applications 
    lo = quantile(imss(:),qlo);
    hi = quantile(imss(:),qhi);
    im_scaled = (im - lo)/(hi-lo);
end
