%% MB Feedback network model
%  Basic model: PNs -> KCs -> MBON-| FBN -> MBIN with a lateral MBN input
%  from an excitatory MBON.

% Common Parameters across KC->MBON weights

p.tmax = 1200;
p.t0 = 500;
p.t1 = 1000;
p.N = 7; % PN, KC, MBON, MBONL, FB, FB2, MBIN
p.tauinv = 1./[10 10 10 10 10 10 10]';
p.V0 = [0, 0, 0, 0, 0, 0, 0]';

q.k = 1.5;
q.yh = 5;
q.rel_strength = [10 10 10 10 10 10 10]';

    
Jpnkc    = 1;
JkcmbonL = 1;
Jmbonfb  = -0.5;
JmbonLfb = 0.5;   % Key parameter in controlling shape of MBIN response.
Jfbfb2   = 0.5;
Jfbmbin  = 0.4;
Jfb2fb   = -0;
Jfb2mbin = -0.5;

tonic_activation = [0 0 0 0 2 2 4]';

% Run simulations across different KC->MBON weights, activating olfaction

Jkcs = 0:0.01:1;
stim = [10, 0, 0, 0, 0, 0, 0]';

mbon_max = zeros(size(Jkcs));
fb_max = zeros(size(Jkcs));
fb2_max = zeros(size(Jkcs));
mbin_max = zeros(size(Jkcs));

mbin_dyn = cell(length(Jkcs),2);
all_dat = cell(length(Jkcs),1);
for jj = 1:length(Jkcs)
    Jkcmbon  = Jkcs(jj);
    
    p.A =   [ [0 0 0 0 0 0 0];...  % PN input
            [Jpnkc 0 0 0 0 0 0];...  % KC input
            [0 Jkcmbon 0 0 0 0 0];...  % MBON input
            [0 JkcmbonL 0 0 0 0 0];...  % MBONL input
            [0 0 Jmbonfb JmbonLfb Jfb2fb 0 0];...  % FB input
            [0 0 0 0 Jfbfb2 0 0];...  % FB2 input
            [0 0 0 0 Jfbmbin Jfb2mbin 0]];    % MBIN input
    
    [r_out, t_out] = logistic_integration_general( stim, p, q, tonic_activation );
    
    mbon_max(jj) = max(r_out(:,3));
    fb_max(jj) = max(r_out(:,5));
    fb2_max(jj) = max(r_out(:,6));
    mbin_max(jj) = max(r_out(:,7));
 
    tss = and( t_out > 850, t_out<900 ); % Steady state sampling window
    mbon_ss(jj) = max(r_out(tss,3));
    fb_ss(jj) = max(r_out(tss,5));
    fb2_ss(jj) = max(r_out(tss,6));
    mbin_ss(jj) = max(r_out(tss,7));
    
    mbin_dyn{jj,1} = t_out;
    mbin_dyn{jj,2} = r_out(:,7);
    all_dat{jj} = r_out;
end


cmap_qual = cbrewer('qual','Dark2',4);
clr1d = 0.9*cbrewer('div','RdBu',length(Jkcs));

% Same thing with steady state values
figure('Color','w'); hold on;
plot(Jkcs,mbon_ss / max(mbon_ss),'Color',cmap_qual(1,:),'LineWidth',2)
plot(Jkcs,fb_ss / max(fb_ss),'Color',cmap_qual(2,:),'LineWidth',2)
plot(Jkcs,fb2_ss / max(fb2_ss),'Color',cmap_qual(4,:),'LineWidth',2)
plot(Jkcs,mbin_ss / max(mbin_ss),'Color',cmap_qual(3,:),'LineWidth',2)
plot([0.1],mbin_ss(0.1==Jkcs) / max(mbin_ss),'Marker','o','MarkerSize',10,'LineWidth',2,'Color',clr1d( 0.1==Jkcs, :) )
plot([0.5],mbin_ss(0.5==Jkcs) / max(mbin_ss),'Marker','o','MarkerSize',10,'LineWidth',2,'Color',clr1d( 0.5==Jkcs, :) )
plot([0.9],mbin_ss(0.9==Jkcs) / max(mbin_ss),'Marker','o','MarkerSize',10,'LineWidth',2,'Color',clr1d( 0.9==Jkcs, :) )
set(gca,'XDir','Reverse','TickDir','out','YTick',[0, 0.5, 1], 'YLim',[0 1] )
axis square
export_fig('steady_state_values_fb_june5.pdf')

mbon_val = zeros(size(Jkcs));
fb_val = zeros(size(Jkcs));
fb2_val = zeros(size(Jkcs));
mbin_val = zeros(size(Jkcs));

figure('Color','w'); hold on;
for ii = 1:length(Jkcs)
    if or( Jkcs(ii) == 0.1, or(Jkcs(ii) == 0.5 , Jkcs(ii) == 0.9 ) )
%        plot(mbin_dyn{ii,1},mbin_dyn{ii,2} / max(max(mbin_max)),'Color',clr1d(ii,:),'LineWidth',2 )
        plot(mbin_dyn{ii,1},mbin_dyn{ii,2},'Color',clr1d(ii,:),'LineWidth',2 )
    end
end
axis square
set(gca,'XLim',[400 1200],'TickDir','out','YLim',[0 5.5]);
export_fig('examples_with_fb2_june5.pdf')

%% Rerun simulation assuming lower/higher baseline FBN states

tonic_activation_low = [0 0 0 0 2 2 4]';
tonic_activation_high = [0 0 0 0 4 2 4]';
Jkcs = 0.9;

stim_mbon = [0, 0, 10, 0, 0, 0, 0]';
p_baseline = p;
p_baseline.A =   [ [0 0 0 0 0 0 0];...  % PN input
        [Jpnkc 0 0 0 0 0 0];...  % KC input
        [0 Jkcmbon 0 0 0 0 0];...  % MBON input
        [0 JkcmbonL 0 0 0 0 0];...  % MBONL input
        [0 0 Jmbonfb JmbonLfb Jfb2fb 0 0];...  % FB input
        [0 0 0 0 Jfbfb2 0 0];...  % FB2 input
        [0 0 0 0 Jfbmbin Jfb2mbin 0]];    % MBIN input

[r_out_low, t_out_low] = logistic_integration_general( stim_mbon, p_baseline, q, tonic_activation_low );
[r_out_high, t_out_high] = logistic_integration_general( stim_mbon, p_baseline, q, tonic_activation_high );

figure('Color','w'); hold on;
plot(t_out_low,r_out_low(:,7),'Color',[0.5, 0.5, 0.5], 'LineWidth',2)
plot(t_out_high,r_out_high(:,7),'Color',[0.1, 0.1, 0.1], 'LineWidth',2)
set(gca,'XLim',[400 1200],'TickDir','out','YLim',[0 4]);
legend({'Low FBN Baseline', 'High FBN Baseline'});
export_fig('baseline_comparison.pdf')

%% Re-do the above at the high baseline values.


Jpnkc    = 1;
JkcmbonL = 1;
Jmbonfb  = -0.5;
JmbonLfb = 0.5;   % Key parameter in controlling shape of MBIN response.
Jfbfb2   = 0.5;
Jfbmbin  = 0.4;
Jfb2fb   = -0;
Jfb2mbin = -0.5;

tonic_activation = tonic_activation_high;

% Run simulations across different KC->MBON weights, activating olfaction

Jkcs = 0:0.01:1;
stim = [10, 0, 0, 0, 0, 0, 0]';

mbon_max = zeros(size(Jkcs));
fb_max = zeros(size(Jkcs));
fb2_max = zeros(size(Jkcs));
mbin_max = zeros(size(Jkcs));

mbin_dyn = cell(length(Jkcs),2);
all_dat = cell(length(Jkcs),1);
for jj = 1:length(Jkcs)
    Jkcmbon  = Jkcs(jj);
    
    p.A =   [ [0 0 0 0 0 0 0];...  % PN input
            [Jpnkc 0 0 0 0 0 0];...  % KC input
            [0 Jkcmbon 0 0 0 0 0];...  % MBON input
            [0 JkcmbonL 0 0 0 0 0];...  % MBONL input
            [0 0 Jmbonfb JmbonLfb Jfb2fb 0 0];...  % FB input
            [0 0 0 0 Jfbfb2 0 0];...  % FB2 input
            [0 0 0 0 Jfbmbin Jfb2mbin 0]];    % MBIN input
    
    [r_out, t_out] = logistic_integration_general( stim, p, q, tonic_activation );
    
    mbon_max(jj) = max(r_out(:,3));
    fb_max(jj) = max(r_out(:,5));
    fb2_max(jj) = max(r_out(:,6));
    mbin_max(jj) = max(r_out(:,7));
 
    tss = and( t_out > 850, t_out<900 ); % Steady state sampling window
    mbon_ss(jj) = max(r_out(tss,3));
    fb_ss(jj) = max(r_out(tss,5));
    fb2_ss(jj) = max(r_out(tss,6));
    mbin_ss(jj) = max(r_out(tss,7));
    
    mbin_dyn{jj,1} = t_out;
    mbin_dyn{jj,2} = r_out(:,7);
    all_dat{jj} = r_out;
end


cmap_qual = cbrewer('qual','Dark2',4);
clr1d = 0.9*cbrewer('div','RdBu',length(Jkcs));

% Same thing with steady state values
figure('Color','w'); hold on;
plot(Jkcs,mbon_ss / max(mbon_ss),'Color',cmap_qual(1,:),'LineWidth',2)
plot(Jkcs,fb_ss / max(fb_ss),'Color',cmap_qual(2,:),'LineWidth',2)
plot(Jkcs,fb2_ss / max(fb2_ss),'Color',cmap_qual(4,:),'LineWidth',2)
plot(Jkcs,mbin_ss / max(mbin_ss),'Color',cmap_qual(3,:),'LineWidth',2)
plot([0.1],mbin_ss(0.1==Jkcs) / max(mbin_ss),'Marker','o','MarkerSize',10,'LineWidth',2,'Color',clr1d( 0.1==Jkcs, :) )
plot([0.64],mbin_ss(0.64==Jkcs) / max(mbin_ss),'Marker','o','MarkerSize',10,'LineWidth',2,'Color',clr1d( 0.5==Jkcs, :) )
plot([0.9],mbin_ss(0.9==Jkcs) / max(mbin_ss),'Marker','o','MarkerSize',10,'LineWidth',2,'Color',clr1d( 0.9==Jkcs, :) )
set(gca,'XDir','Reverse','TickDir','out','YTick',[0, 0.5, 1], 'YLim',[0 1] )
axis square
export_fig('high_baseline_ss_values_june5.pdf')
