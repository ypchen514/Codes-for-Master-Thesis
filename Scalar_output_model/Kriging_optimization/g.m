function y = g(alpha,p1,p2,p3,p4,matrix)
for i = 1:size(matrix,1)
    R = matrix(i,:);
    y(i) = alpha+exp(p1*((R(1)-1)^2+R(2)^2)+p2*(R(3)^2+R(4)^2))+exp(p3*(R(1)^2+(R(2)-1)^2)+p4*(R(3)^2+R(4)^2));
end

end