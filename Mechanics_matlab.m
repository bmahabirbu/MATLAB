%Github version
%%X and Y vectors are row vectors not column vectors
load('TrussDesign.mat','C','Sx','Sy','X','Y','L');

%%code for formating the A and T matrix
[jointNum, memberNum] = size(C); %%connection matrix C. This matrix has j rows and m columns
Cx = zeros(jointNum, memberNum); 
Cy = zeros(jointNum, memberNum);
[jointLocation, memberLocation] = find(C == 1);
 
%%Make Cx ie find where C = 1 and do sum of Fx
for i = 1:(length(jointLocation)/2)
   vec = find(C(:,i) == 1);
   Cx(jointLocation(2*i-1), memberLocation(2*i-1)) = (X(vec(2))-X((vec(1))))/sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
   Cx(jointLocation(2*i), memberLocation(2*i)) = (X(vec(1)) -X((vec(2))))/sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
end

%%Make Cy ie find where C = 1 and do sum Fy
for i = 1:(length(jointLocation)/2)
   vec = find(C(:,i) == 1);
   Cy(jointLocation(2*i-1), memberLocation(2*i-1)) = (Y(vec(2))-Y((vec(1))))/sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
   Cy(jointLocation(2*i), memberLocation(2*i)) = (Y(vec(1)) -Y((vec(2))))/sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2)); 
end

%%add matrix and create A
A = [Cx Sx;Cy Sy];
 
%%code for finding the unknown value matrix T
T = A\L;

%%code for printing unknown force in Newtons
fprintf('EK301, Section A3, Group: Harin Lee, Landon Kushimi, Brian Mahabir \n');
fprintf('Load: %d\n',L(L ~= 0));
%disp(memberNum); %%debug
fprintf('Member forces in Newtons:\n');
for i = 1:memberNum
    if round(T(i),3) == 0
        fprintf('m%d: Zero force member\n',i);
    elseif round(T(i),3) > 0
        fprintf('m%d: %.3f (T)\n',i,abs(T(i)));
    else 
         fprintf('m%d: %.3f (C)\n',i,abs(T(i)));
    end
end

%%Code for printing reaction forces
fprintf('Reaction forces in Newtons:\n');
fprintf('Sx1: %.3f \n',round(T(end-2),2,'significant'));
fprintf('Sy1: %.3f \n',round(T(end-1),2,'significant'));
fprintf('Sy2: %.3f \n',round(T(end),2,'significant'));
totalLength = 0;

%%code for cost
for i = 1:memberNum
            vec = find(C(:,i) == 1);
            totalLength = totalLength + sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
end
cost = 10 * jointNum + totalLength;
fprintf('Cost of truss: $%.2f\n', cost);

%%code for finding theoretical max value
%%create new variables as to not change values old ones (for debugging purposes)
nL = L;
nSx = Sx;
nSy = Sy;

Load = find(nL(:,1) > 0);
nL(Load,1) = 1; %%add 1 Newton value to load node

memberlength = zeros(memberNum,1); %%create a matrix of member lengths

%%code for formating the A and T matrix
[jointNum, memberNum] = size(C); %%connection matrix C. This matrix has j rows and m columns
nCx = zeros(jointNum, memberNum); 
nCy = zeros(jointNum, memberNum);
[jointLocation, memberLocation] = find(C == 1);
 
%%Make nCx ie find where C = 1 and do sum of Fx
for i = 1:(length(jointLocation)/2)
   vec = find(C(:,i) == 1);
   memlen = sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2)); %%get length bewtween joints
   memberlength(i,1) = memlen; %%put lengths in array 
   nCx(jointLocation(2*i-1), memberLocation(2*i-1)) = (X(vec(2))-X((vec(1))))/memlen;
   nCx(jointLocation(2*i), memberLocation(2*i)) = (X(vec(1)) -X((vec(2))))/memlen;
end

%%Make nCy ie find where C = 1 and do sum Fy
for i = 1:(length(jointLocation)/2)
   vec = find(C(:,i) == 1);
   memlen = sqrt((((X(vec(2))-X(vec(1))))^2) + ((Y(vec(2))-Y(vec(1)))^2));
   nCy(jointLocation(2*i-1), memberLocation(2*i-1)) = (Y(vec(2))-Y((vec(1))))/memlen;
   nCy(jointLocation(2*i), memberLocation(2*i)) = (Y(vec(1)) -Y((vec(2))))/memlen; 
end

%%add matrix and create A
nA = [nCx nSx; nCy nSy];
 
%%code for finding the unknown value matrix T
nT = nA\nL;

%%Uncertainty code using Euler fit (also could use Empirical still havent
%%decided
Upper_Confidence = 1340.9./(memberlength.^2) + 20.6./(memberlength.^2); 
Lower_Confidence = 1340.9./(memberlength.^2) - 20.6./(memberlength.^2);
Buckling_load = (Upper_Confidence+Lower_Confidence)/2;

%%resize nT to not include support forces
[Lent, Widt] = size(nT);
nT(Lent-2:Lent) = []; %%delete last 3 forces which are support 

%%resize solution matrix to only have compression members and resize (also
%%resize memberlength matrix accordingly
compression = find(nT(:,1) > 0);
[Len, Width] = size(compression);
for i = 1 : Len
    nT(compression(i,1)) = 0;
    memberlength(compression(i,1)) = 0;
end

%%find SR using the buckling force
Sr = abs(nT./Buckling_load);

Sr_max = max(Sr(:));
Failure = 1/Sr_max;
fprintf('max theoretical load is: %f N\n',Failure);

%fprintf('Theoretical max load/cost ratio in N/$: $%.5f\n', maxRatio);
