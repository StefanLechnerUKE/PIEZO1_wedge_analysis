
function [result] = PlotTrimerAnalysisResult(result)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot analysis results
if size(result.Trimers,1) > 0
screenSize = get(0, 'ScreenSize');
figure('Position',[0 screenSize(4) 900 300]);
tiledlayout(1,6);

%% plot mean interblade
nexttile
    result.InterBlades = cat(2, result.InterBlades, mean(result.InterBlades(:,1:3),2));
    bar(mean(result.InterBlades(:,4)),'FaceColor',result.color);
    hold on;
    swarmchart(ones(size(result.InterBlades,1)),result.InterBlades(:,4),'k','filled');
    ylabel('interblade distance(nm)');

%% plot NN vs trimer interblade
    Meaninter = mean(result.InterBlades(:,1:3),2);
    nexttile(2,[1,2])
    scatter(result.InterBladesWithNN(:,4),Meaninter,100,result.color, 'filled');
    xlim([0 500]); ylim([0 50]);
    title('interblade vs. nearest neighbor'); xlabel('distance to nearest neighbor (nm)'); ylabel('interblade distance (nm)');

%% plot interblade angle histogram
    result.AnglesSingleColumn = cat(1, result.Angles(:,1), result.Angles(:,2),result.Angles(:,3));
    nexttile
    histogram(result.AnglesSingleColumn, [0:10:120],'FaceColor',result.color);
    title('Interblade angle distribution'); xlabel('interblade angle (°)'); ylabel('counts');

%% plot superparticle
[result.Superparticle] = PIEZO1Superparticle(result); % second parameter: choose 1 align or no align
nexttile(5,[1,2])  
hold on;
    scatter(result.Superparticle(:,1), result.Superparticle(:,2),200, "o", "filled", 'MarkerFaceColor', result.color, 'MarkerFaceAlpha',[0.3]);
    scatter(result.Superparticle(:,3), result.Superparticle(:,4),200, "o", "filled", 'MarkerFaceColor', result.color, 'MarkerFaceAlpha',[0.3]);
    scatter(result.Superparticle(:,5), result.Superparticle(:,6),200, "o", "filled", 'MarkerFaceColor', result.color, 'MarkerFaceAlpha',[0.3]);
    grid on; axis equal; xlim([-30 30]); ylim([-30 30]);
end
