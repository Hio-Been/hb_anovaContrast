function [pval, ci] = hb_anovaContrast(dataset, cont_coeffs, alpha)
%% INPUT data type
% 
% dataset <- cell([ 1, nCondition ]); % each cell contains column vector
% cont_coeffs <- [ 1, -1, 0, 0 ...]; % should be length at nCondition
% 
% written by Hio-Been Han, for GOSHIMTONG 2016-1, YONSEI PSYCHOLOGY
% hiobeen@yonsei.ac.kr, 20160420

if nargin < 3
    alpha = .05;
end

%% START
disp(['.']);
disp(['.']);
disp(['-----hb_anovaContrast.m']);

if ~iscell(dataset) % Same N for each condition
    [ nSbj_perCond, nCond ] = size(dataset);
    dat = cell([1, nCond]);
    for i = 1:nCond
        dat{1,i} = dataset(:,i);
    end
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
for i = 1:nCond
    Ns(i) = length(dat{1,i});
    Ms(i) = mean(dat{1,i});
    Ss(i) = std(dat{1,i},0);
    disp([ 'cond ' num2str(i) '_N: ' num2str(Ns(i)) ', M : ' num2str(Ms(i)) ', S : ' num2str(Ss(i))]);
end

%% Sp
Sp_square =sum( (Ns-1) .* (Ss.^2)) / sum((Ns-1));
Sp = sqrt(Sp_square);
disp(['Sp : ' num2str(Sp) ', Sp_square : ' num2str(Sp_square) ]);
disp(['-----']);

%% Contrast test
disp(['Weights : [ ' num2str(cont_coeffs) ' ]' ]);
md = sum(Ms .* cont_coeffs);
se = sqrt( sum( ((cont_coeffs.^2)*Sp_square) ./ Ns ) );
df = sum( Ns-1 );
t_crit = abs( tinv( alpha*.5, df) ); 
disp(['t_crit (original) : ' num2str(t_crit) ]);
t_crit_round = round(t_crit*100)/100;
disp([' rounded t_crit : ' num2str(t_crit_round)]); % for GOSHIMTONG policy
disp(['-----']);
disp(['-----RESULT----']);
disp(['-----']);

%% Disp Confidence Interval
ci =md + ([-1 1]*(t_crit_round*se));
disp([' MD : ' num2str(md)  ', SE : ' num2str(se) ', t_crit: ' num2str(t_crit_round) ', df : ' num2str(df)])
disp([' CI : [ ' num2str(ci(1)) ', ' num2str(ci(2)) ' ]']);

t_stat = md / se;
pval = 1-tcdf(t_stat, df);
if pval < alpha
    result_txt = 'Reject H0. diff is real';
else
    result_txt = 'Null result'
end
disp([ result_txt ', t(' num2str(df) ') = ' num2str(t_stat) ', p = ' num2str(pval) ]);

return




