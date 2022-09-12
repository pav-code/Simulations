function b = backProjection2(g2,theta)
b = zeros(774);
for x=-386:387
    for y=-386:387
        b(x+387,y+387) = g2(548+round(x*cosd(theta)+y*sind(theta)), theta);
    end
end
end