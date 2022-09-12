function b = backProjection(g,theta)
b = zeros(100);
for x=-49:50
    for y=-49:50
        b(x+50,y+50) = g(72+round(x*cosd(theta)+y*sind(theta)), theta);
    end
end
end

