function confusion_matrix(act,det)

[mat,order] = confusionmat(act,det);
%mat = mat/length(act);
k=max(order);             %k is the number of class

imagesc(mat); %# Create a colored plot of the matrix values
colormap(flipud(gray));  %# Change the colormap to gray (so higher values are

%#black and lower values are white)
title('Confusion matrix'); 
textStrings = num2str(mat(:),'%0.02f');       %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding

%# Create x and y coordinates for the strings 
[x,y] = meshgrid(1:k);  
hStrings=text(x(:),y(:),textStrings(:),'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(mat(:) > midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color

set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors；


%
set(gca,'XTick',0:k,...                                    
        'YTick',0:k,...                                    %
        'TickLength',[0 0]);
xlabel('Predicted');
ylabel('Actual');

overall = 0;
rowMulcol = 0;
for i =1:k
    overall = overall + mat(i,i);
    rowMulcol = rowMulcol+sum(mat(i,:))*sum(mat(:,i));
    fprintf('第%d类：\n', i);
    fprintf('制图精度=%8.5f',mat(i,i)/sum(mat(i,:)));
    fprintf('  漏分误差=%8.5f\n',1-mat(i,i)/sum(mat(i,:)));
    fprintf('用户精度=%8.5f',mat(i,i)/sum(mat(:,i)));
    fprintf('  错分误差=%8.5f \n',1-mat(i,i)/sum(mat(:,i)));
    fprintf('\n');
end
overall = overall/length(act);
rowMulcol = rowMulcol/(length(act).^2);
kappa = (overall-rowMulcol)/(1-rowMulcol);
fprintf('总体精度=%8.5f\n',overall);
fprintf('Kappa=%8.5f\n',kappa);


