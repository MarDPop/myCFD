Table  = dlmread('RocketCFD\meshFaces.csv',',');
figure
hold on
for i = 1:length(Table(:,1))
    plot([Table(i,1) Table(i,3)],[Table(i,2) Table(i,4)],'b');
end
axis equal