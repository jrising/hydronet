function flows = model(elevation, glaciers, mu_s, mu_i, umask, ...
                       precips, prows, u2p, temps, trows, u2t)
% TO FIX: can't have convolution every day-- totally different
% path.  Collect all, then do once.
[rows, cols] = size(elevation);
[pall, pcols] = size(precips);
[tall, tcols] = size(temps);

snow = zeros(rows, cols);
flows = zeros(1, min(pall / prows, tall / trows));

dayflows = zeros(1, floor(max(max(umask)) + 1));

for tt = 1:length(flows)
  tempstt = scale(ldeott(temps, tt, trows), u2t);
  precipstt = scale(ldeott(precips, tt, prows), u2p);

  melt_i = max(glacier, mu_i(elevation) .* tempstt .* (glaciers > 0) .* (tempstt > 0) .* (snow <= 0));
  glacier = glacier - melt_i;
  melt_s = max(snow, mu_s(elevation) .* tempstt .* (tempstt > 0) .* (snow > 0));
  snow = snow + precipstt .* (tempstt < 0) - melt_s;
  
  flow = (precips + melt_s + melt_i) .* (tempstt > 0);
  for dd = 1:length(dayflows)
    dayflows(dd) = dayflows(dd) + sum(sum(flow(umask >= (dd - 1) && umask < dd)));
  end
  
  dayflows = conv([.5 .5], dayflows);
  flows(tt) = dayflows(1);
  dayflows = dayflows(2:end);
end
  