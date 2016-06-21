function colorInterval = findColorInterval(image1,image2)

colMax = max(max(max(image1)),max(max(image2)));
colMin = min(min(min(image1)),min(min(image2)));

colorInterval = [colMin colMax];

end