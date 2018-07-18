function a=mslide(S,loc,num);
%mslide(S,loc,num); matrix slide: deletes k rows after (including) loc, and returns the trimmed matrix

part1 = S(:,1:loc-1);

part2 = S(:,loc+num:size(S,2));
a = [part1 part2];
