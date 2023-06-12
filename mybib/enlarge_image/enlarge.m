function I1 = enLarge(I,LeftUpPoint, RightBottomPoint, Enlargement_Factor, LineWidth)
%I = imread('cameraman.tif');
%I = double(I);
%I = I/max(I(:));
%I = Re_tensor{10};
% I1 = ShowEnlargedRectangle(I,[100 119], [166 186] ,1.6, 2);
I1 = ShowEnlargedRectangle(I,LeftUpPoint, RightBottomPoint, Enlargement_Factor, LineWidth);
% [63 134], [96 165] [129 106], [170 141] [24 58], [86 132]
% [1 19], [32 60] [129 106], [170 141] [24 58], [86 132]
% 55*55
% [77 163], [138 221] , 1.9 house; [29 59], [84 114] 2 tulips;
% [1 20], [39 59] , 2.9 barbara; [88 120], [163 193] , 1.5 pepper;
% [94 119], [166 186] ,1.6 lena
% [131 143], [168 181] ,2.9 toy; beads [173 160], [236 224] ,1.7 [120 195], [166 243] ,2.3
% [78 102], [147 172] ,1.5 cloth;
% figure;
imshow(I1,'border','tight','initialmagnification','fit');
set(gcf,'Units','centimeters');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition([4,4])],'PaperSize',screenposition([4,4]));
%print('test.pdf','-dpdf');
