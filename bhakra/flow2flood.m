function [vols, iis] = flow2flood(flows, year0)

years = floor(length(flows) / 365);
days = daysperyear(year0:year0+years);
indexes = cumsum([1 days]);

% 1. Determine 2 year peak flow
%  - find peak flow in each year
maxes = [];
for ii = 1:length(indexes) - 1
  if indexes(ii+1) > length(flows)
    break;
  end
  maxes = [maxes max(flows(indexes(ii):indexes(ii+1)))];
end

%  - find the median: this is the 2 year peak flow
threshold = median(maxes);

% 2. Determine cummulative volume over threshold
vols = [];
iis = [];
flooding = 0;
for ii = 1:length(flows)
  if flooding
    if flows(ii) >= threshold
      flooding = flooding + flows(ii) - threshold;
    else
      vols = [vols flooding];
      flooding = 0;
    end
  else
    if flows(ii) >= threshold
      flooding = flows(ii) - threshold;
      iis = [iis ii];
    end
  end
end

figure
plot(flows);
hold on
plot([1 length(flows)], [threshold threshold])
plot(iis, threshold * ones(1, length(iis)), '.r');
yearaxis(gca, 10, 1962, 2008);
axis tight