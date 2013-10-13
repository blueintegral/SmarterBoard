function show_digit(img)

img_copy = img .- min(img);
img_copy = (256/max(img_copy)).*img_copy;

colormap(1 .- gray);
image(reshape(img_copy, 16, 16)');


jmee@umiche