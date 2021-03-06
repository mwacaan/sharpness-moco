FWHM=0.7*2.355; % conversion sigma to FWHM, given imaging resolution of 0.7 mm


%% calculate some pre-stats
CRthresh=pars.mrelcb;
% first corr, then uncorr
for ii=1:length(AR)
  if ~isempty(AR{ii})
    % compute CRLB
    cbcorr=2*AR{ii}.errcorr./AR{ii}.sigmacorr;
    cbuncorr=2*AR{ii}.erruncorr./AR{ii}.sigmauncorr;
    
    medcb(ii,:)=[nanmedian(cbuncorr) nanmedian(cbcorr)];

    % update aggregate test score (binary)
    tmp=abs(AR{ii}.aggregates{1}); 
    tmp4=tmp>=4;
    tmp=mod(tmp,4);
    test=abs(cbcorr)>CRthresh;
    agg1=mod(tmp,2)-2*test-4*tmp4;
    tmp=abs(AR{ii}.aggregates{2}); 
    tmp4=tmp>=4;
    tmp=mod(tmp,4);
    test=abs(cbuncorr)>CRthresh;
    agg2=mod(tmp,2)-2*test-4*tmp4;
    validboth{ii}=and(~agg1,~agg2);
    
    sc = AR{ii}.sigmacorr;
    ec = AR{ii}.errcorr;
    su = AR{ii}.sigmauncorr;
    eu = AR{ii}.erruncorr;
    unc = sqrt(abs(1./(su.^2-sc.^2).*(su.^2.*eu.^2+sc.^2.*ec.^2)));

    sigma_diff_abs = validboth{ii}.*sqrt(abs(su.^2-sc.^2));
    signed_indication = 2.*((su.^2-sc.^2)>0)-1; % 1 if better, -1 if worse
    sigma_diff_signed{ii} = signed_indication.*sigma_diff_abs;

    better_signed{ii} = validboth{ii}.*sigma_diff_signed{ii}.*((su-eu)>(sc+ec));
    worse_signed{ii} = validboth{ii}.*sigma_diff_signed{ii}.*((sc-ec)>(su+eu));
    uncertainty{ii} = validboth{ii}.*unc;
  end
end

nROI=length(mot_mean)/nsubj;
% relative nr of valid clusters
relvalid=cellfun(@sum,validboth)./cellfun(@length,validboth);
relvalid=reshape(relvalid,nROI,[]);

% number of clusters
ln=ones(length(AR),1); ln(:)=nan;
for ii=1:length(AR)
  try
  ln(ii)=length(AR{ii}.sigmacorr);
  end
end
ln=reshape(ln,nROI,[]);


allsigma = []; allsigmacorr=[]; allsigmauncorr=[];
allmotion = []; allu=[];
subj_name = cell(1,length(doSubj));
clear csigma cu
plotIDs=[];
i=iSubj;
  subj = find(~cellfun(@isempty,strfind(fields,FileID.uIDs{i})));
  plot_sigma = []; plot_sigma_corr=[]; plot_sigma_uncorr=[];
  plot_motion = []; plot_u=[];
  for j = 1:nROI
      inds = j;
      fieldinds = subj(inds);
      sigm_plot = []; sigm_mean_plot=[]; 
      sigm_corr=[]; sigm_uncorr=[];
      mot_plot = []; mot_mean_plot=[];
      k = 1;
      fieldind = fieldinds(k);
      sigm_data = better_signed{fieldind}+worse_signed{fieldind};
      signif = better_signed{fieldind}|worse_signed{fieldind};
      tmp=sigma_diff_signed{fieldind}(validboth{fieldind}>0);
      u=uncertainty{fieldind}(validboth{fieldind}>0);
      try
        csigma{i}=cat(2,csigma{i},FWHM*tmp);
        cu{i}=cat(2,cu{i},FWHM*u);
      catch
        csigma{i}=FWHM*tmp;
        cu{i}=FWHM*u;
      end
      sigm_mean_plot(k)=FWHM*median(tmp);
      uplot(k)=FWHM*nanmedian(u);
      tmp=AR{fieldind}.sigmacorr(validboth{fieldind}>0);       % sigma after
      sigm_corr(k)=FWHM*median(tmp);
      tmp=AR{fieldind}.sigmauncorr(validboth{fieldind}>0);       % sigma after
      sigm_uncorr(k)=FWHM*median(tmp);
      mot_mean_plot(k) = (mot_mean{fieldind});
      mot_add = ones(sum(signif),1)*mot_mean{fieldind}; 
      mot_plot = cat(2,mot_plot,mot_add');
      plot_sigma = cat(2,plot_sigma,sigm_mean_plot);
      plot_u=cat(2,plot_u,uplot);
      plot_motion = cat(2,plot_motion,mot_mean_plot);
      plot_sigma_corr = cat(2,plot_sigma_corr,sigm_corr);
      plot_sigma_uncorr = cat(2,plot_sigma_uncorr,sigm_uncorr);
  end
  allsigma = cat(2,allsigma,plot_sigma);
  allu=cat(2,allu,plot_u);
  allsigmacorr = cat(2,allsigmacorr,plot_sigma_corr);
  allsigmauncorr = cat(2,allsigmauncorr,plot_sigma_uncorr);
  allmotion = cat(2,allmotion,plot_motion);
  
csigma=csigma(doSubj); cu=cu(doSubj);

%%

disp('ROIs:')
disp(fields')

disp('FWHM uncorrected:')
disp(num2str(allsigmauncorr))
disp('FWHM corrected:')
disp(num2str(allsigmacorr))
disp('FWHM difference:')
disp(num2str(allsigma))

disp('Number of valid clusters:')
disp(num2str(ln'))

disp('Relative nr of valid clusters:')
disp(num2str(relvalid'))

