function [distance] = dist2(coord1, coord2)

distance = sqrt((coord1(1)-coord2(1))^2 + (coord1(2)-coord2(2))^2);
