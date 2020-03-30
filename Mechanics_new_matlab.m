load('TrussDesign.mat','C','Sx','Sy','X','Y','L');
[jointNum, memberNum] = size(C);
Cx = zeros(jointNum, memberNum);
Cy = zeros(jointNum, memberNum);
[jointLocation, memberLocation] = find(C == 1);
for i = 1:(length(jointLocation)/2)
   vec = find(C(:,i) == 1);
   Cx(jointLocation(2*i-1), memberLocation(2*i-1)) = (X(vec(2))-X((vec(1))))/sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
   Cx(jointLocation(2*i), memberLocation(2*i)) = (X(vec(1)) -X((vec(2))))/sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
end
for i = 1:(length(jointLocation)/2)
   vec = find(C(:,i) == 1);
   Cy(jointLocation(2*i-1), memberLocation(2*i-1)) = (Y(vec(2))-Y((vec(1))))/sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
   Cy(jointLocation(2*i), memberLocation(2*i)) = (Y(vec(1)) -Y((vec(2))))/sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2)); 
end
A = [Cx Sx;Cy Sy];
 
T = A\L;
fprintf('EK301, Section A3, Group: Harin Lee, Landon Kushimi, Brian Mahabir \n');
fprintf('Load: %d\n',L(find(L ~= 0)));
fprintf('Member forces in Newtons:\n');
for i = 1:memberNum
    if T(i) > 0
        fprintf('m%d: %.3f (T)\n',i,abs(T(i)));
    elseif T(i) == 0
        fprintf('m%d: Zero force member\n',i);
    else
        fprintf('m%d: %.3f (C)\n',i,abs(T(i)));
    end
end
fprintf('Reaction forces in Newtons:\n');
fprintf('Sx1: %.3f \n',T(end-2));
fprintf('Sy1: %.3f \n',T(end-1));
fprintf('Sy2: %.3f \n',T(end));
totalLength = 0;
for i = 1:memberNum
            vec = find(C(:,i) == 1);
            totalLength = totalLength + sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
end
cost = 10 * jointNum + totalLength;
fprintf('Cost of truss: $%.2f\n', cost);
maxRatio = (L(find(L ~= 0)))/cost;
fprintf('Theoretical max load/cost ratio in N/$: $%.5f\n', maxRatio);
