a=simulation_bending_protein(8);
local8_l1 = a(1,:);

b=simulation_bending_protein(9);
local9_l1=b(1,:);

c=simulation_bending_protein(10);
local10_l1 = c(1,:);

d=simulation_bending_protein(11);
local11_l1=d(1,:);

e=simulation_bending_protein(12);
local12_l1 = e(1,:);

f=simulation_bending_protein(13);
local13_l1 = f(1,:);

g=simulation_bending_protein(14);
local14_l1 = g(1,:);

xlim = length(local8_l1);

local8_l1N = [];
local9_l1N = [];
local10_l1N = [];
local11_l1N = [];
local12_l1N = [];
local13_l1N = [];
local14_l1N = [];

for i = 1:xlim
    local8_l1N(i) = local8_l1(xlim - i + 1);
    local9_l1N(i) = local9_l1(xlim - i + 1);
    local10_l1N(i) = local10_l1(xlim - i + 1);
    local11_l1N(i) = local11_l1(xlim - i + 1);
    local12_l1N(i) = local12_l1(xlim - i + 1);
    local13_l1N(i) = local13_l1(xlim - i + 1);
    local14_l1N(i) = local14_l1(xlim - i + 1);
end

plot(1:xlim,local8_l1N,'b','LineWidth',2)
hold on
plot(1:xlim,local9_l1N,'g','LineWidth',2)
hold on
plot(1:xlim,local10_l1N,'r','LineWidth',2)
hold on
plot(1:xlim,local11_l1N,'m','LineWidth',2)
hold on
plot(1:xlim,local12_l1N,'k','LineWidth',2)
hold on
plot(1:xlim,local13_l1N,'c','LineWidth',2)
hold on
plot(1:xlim,local14_l1N,'y','LineWidth',2)
hold off

Ang=char(197);
%title('')
set(gca,'FontSize',20)
xlabel('Protein frame number (after renumbering)')
ylabel('Total inconsistency (TIs)')
legend(strcat('Threshold 8',Ang),strcat('Threshold 9',Ang),strcat('Threshold 10',Ang),strcat('Threshold 11',Ang),strcat('Threshold 12',Ang),strcat('Threshold 13',Ang),strcat('Threshold 14',Ang),'Location','northwest')