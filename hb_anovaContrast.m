function [pval, ci, stats] = hb_anovaContrast(dataset, cont_coeffs, alpha, dispOption, onetailOption )
%% INPUT data type
% 
% dataset <- cell([ 1, nCondition ]); % each cell contains column vector
% cont_coeffs <- [ 1, -1, 0, 0 ...]; % should be length at nCondition
% 
% written by Hio-Been Han, for GOSHIMTONG 2016-1, YONSEI PSYCHOLOGY
% hiobeen@yonsei.ac.kr, 20160420 initial
% hiobeen.han@kaist.ac.kr, 20170914 adding some options
% 
if nargin < 3
    alpha = .05;
end
if nargin < 4
    dispOption = true;
end
if nargin < 5
    onetailOption = 0;
end

%% START
% disp(['..']);
% disp(['...']);
% disp(['-----hb_anovaContrast.m']);
% 
if ~iscell(dataset) % Same N for each condition
    [ nSbj_perCond, nCond ] = size(dataset);
    dat = cell([1, nCond]);
    for i = 1:nCond
        dat{1,i} = dataset(:,i);
    end
    % anova1(dat)
else % Generalization for Different N in each condition
    dat = dataset;
    nCond = size(dat, 2);
end
if ~(length(cont_coeffs) == nCond)
    error('Weight LENGTH is incorrect');
end
if ~(sum(cont_coeffs) == 0)
    error('Weight SUM is non-zero');
end

Ns = nan([1, nCond]);
Ms = nan([1, nCond]);
Ss = nan([1, nCond]);
dat_for_anovaN =[];
for i = 1:nCond
    Ns(i) = length(dat{1,i});
    Ms(i) = nanmean(dat{1,i});
    Ss(i) = nanstd(dat{1,i},0);
    if dispOption; disp([ 'cond ' num2str(i) '_N: ' num2str(Ns(i)) ', M : ' num2str(Ms(i)) ', S : ' num2str(Ss(i))]); end
    for i2 = 1:Ns(i)
        dat_for_anovaN = [dat_for_anovaN; [dat{1,i}(i2) ,i ]];
    end
end

stats = [];
stats.Ns = Ns;
stats.Ms = Ms;
stats.Ss = Ss;

%% Sp (pooled variance)
Sp_square =sum( (Ns-1) .* (Ss.^2)) / sum((Ns-1));
Sp = sqrt(Sp_square);
if dispOption
disp(['Sp : ' num2str(Sp) ', Sp_square : ' num2str(Sp_square) ]);
disp(['-----']); end

%% Contrast test
md = sum(Ms .* cont_coeffs);
se = sqrt( sum( ((cont_coeffs.^2)*Sp_square) ./ Ns ) );
df = sum( Ns-1 );
t_crit = abs( tinv( alpha*.5, df) ); 
t_crit_round = round(t_crit*100)/100;
%% Disp Confidence Interval
ci =md + ([-1 1]*(t_crit_round*se));
if dispOption
disp(['Weights : [ ' num2str(cont_coeffs) ' ]' ]);
disp(['t_crit (original) : ' num2str(t_crit) ]);
disp([' rounded t_crit : ' num2str(t_crit_round)]); % for the GOSHIMTONG policy
disp(['-----']);
disp(['-----RESULT----']);
disp([' MD : ' num2str(md)  ', SE : ' num2str(se) ', t_crit: ' num2str(t_crit_round) ', df : ' num2str(df)])
disp([' CI : [ ' num2str(ci(1)) ', ' num2str(ci(2)) ' ]']);
end

t_stat = md / se;
pval = 1-tcdf(abs(t_stat), df);
if dispOption; 
if pval < alpha
    result_txt = '*** Reject H0';
else
    result_txt = 'Null. Go with H0'
end
end
if onetailOption
    pval = pval*.5;
end
if dispOption; 
    disp([ result_txt ', t(' num2str(df) ') = ' num2str(t_stat) ', p = ' num2str(pval) ]);
end
stats.t_stat = t_stat;
stats.df = df;

%% Optional : Multiple comparisons
% [p,table,stats,terms]=anovan(dat_for_anovaN(:,1), dat_for_anovaN(:,2), 'model', 'interaction', 'display', 'on');
% [c,m,h,nms] = multcompare(stats); open c;

return


