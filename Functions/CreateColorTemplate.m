function [colorspec1,colorspec2]=CreateColorTemplate(n);

colors=distinguishable_colors(n);
colors2=[colors,ones(n,1).*0.5];

for i=1:n
    colorspec1(i).spec=colors(i,:);
    colorspec2(i).spec=colors2(i,:);
end

end
