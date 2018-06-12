clear all
close all
parentDir = '~/Bethany/paclab';
% parentDir = '~/code/pac/paclab';
addpath(genpath(parentDir));
results = load(sprintf('%s/Subject_folders/1_block1/threshold.mat',parentDir));
block2 = load(sprintf('%s/Subject_folders/1_block2/threshold.mat',parentDir));

stimulusReversal = results.stimulusReversal;
nReverse = results.nReverse;
%%

sumReversal(1) = sum(stimulusReversal(1,4:nReverse(1)));
sumReversal(2) = sum(stimulusReversal(2,4:nReverse(2)));

stair1mean = sumReversal(1)/(nReverse(1)-3);
stair2mean = sumReversal(2)/(nReverse(2)-3);

StandardDev1 = std(stimulusReversal(1,4:nReverse(1)));
StandardDev2 = std(stimulusReversal(2,4:nReverse(2)));

plot(stimulusReversal(1,1:nReverse(1)));hold on;
plot(stimulusReversal(2,1:nReverse(2)));hold on;

%% Plotting the staircase
%
% It's useful to look at how the intensity (coherence) values changed from
% trial to trial.  This can be seen by simply plotting the values in
% 'results.intensity' with matlab's 'stairs' function.  We'll plot the
% intensity values on a log axis so that they're equally spaced.

figure(1)
clf
stairs(log(results.r1));

%%
% Let's get fancy and plot green and red symbols where trials had correct
% and incorrect responses respectively:

correctTrials = results.response==1;
hold on
plot(find(correctTrials),log(results.intensity(correctTrials)),'ko','MarkerFaceColor','g');

incorrectTrials = results.response==0;
plot(find(incorrectTrials),log(results.intensity(incorrectTrials)),'ko','MarkerFaceColor','r');
set(gca,'YTick',log(2.^[-4:0]))
logy2raw;

xlabel('Trial Number')
ylabel('Coherence');


%%
% You can see from my example (running myself sitting in my office) that
% the staircase is hovering around a coherence level between 0.1 and 0.2.
%
% Let's plot the psychometric function from the same data using similar
% code as before.  The one change we'll make is to have the size of each
% symbol grow with the number of trials presented at the corresponding
% stimulus intensity.  This lets us get a feeling for the 'weight' of each
% data point.  (This wasn't necessary for constant stimuli.  Why not?)

% intensities = unique(results.intensity);
intensities = unique(abs(results.r1)); % all the dif intensity values, or ori
intensities = intensities(2:end);
acc = results.acc;

%
% Then we'll loop through these intensities calculating the proportion of
% times that 'response' is equal to 1:

nCorrect = zeros(1,length(intensities));
nTrials = zeros(1,length(intensities));

for i = 1:length(intensities)
    id = abs(results.r1) == intensities(i);
    nTrials(i) = sum(id);
    nCorrect(i) = sum(results.acc(id));
end

pCorrect = nCorrect./nTrials;

%%

figure(2)
clf
hold on
 plot(log(intensities),100*pCorrect,'-','MarkerFaceColor','b');
 %loop through each intensity so each data point can have it's own size.
for i=1:length(intensities)
    sz = nTrials(i)+2;
    plot(log(intensities(i)),100*pCorrect(i),'ko-','MarkerFaceColor','b','MarkerSize',sz);
end

set(gca,'XTick',log(intensities));
logx2raw;
set(gca,'YLim',[0,100]);
xlabel('Coherence');
ylabel('Percent Correct');

%% NON LOG

figure(2)
clf
hold on
plot(intensities,100*pCorrect,'-','MarkerFaceColor','b');

%loop through each intensity so each data point can have it's own size.
for i=1:length(intensities)
    sz = nTrials(i)+2;
    plot(intensities(i),100*pCorrect(i),'ko-','MarkerFaceColor','b','MarkerSize',sz);
end

set(gca,'XTick',intensities);
set(gca,'YLim',[0,100]);
xlabel('Coherence');
ylabel('Percent Correct');