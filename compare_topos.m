% this function checks if two topographies are similar by means of a 
% permutation test of normalized (GFP) differences. It is based on a
% script written by Hadiseh Nowparast Rostami. I have re-written the script
% into a vectorized version for performance reasons. Also I improved the
% creation of the permutation distribution. Lastly, I have added plotting 
% functionality and encapsulated everything into a function.
% -------------------------------------------------------------------------
% @author: Florian Niefind
% @date: 2016-01-05
% @contact: florian.niefind@posteo.de
% @note: needs eeglab active for topoplot function
% @note: some of the comments assume that data are average-referenced, but
% the code works with any reference!

function [same, p_value, real_diff, f, x]  = compare_topos(topos_1, topos_2, runs, plot, chanlocs, verbose)

    %format: channels * subjects!!
    if verbose
        fprintf('\nNumber of Channels: %i', size(topos_1,1))
        fprintf('\nNumber of Subjects: %i', size(topos_1,2))
    end
    
    %% preparations
    %todo: make these editeable as parameters someday
    maplim = [-2 2];
    
    % topoplot radius
    plotrad = .70;
    headrad = .69;
    intrad  =  1;

    %% compute real differences
    
    %grand average
    topo_1 = mean(topos_1,2);
    topo_2 = mean(topos_2,2);

    %GFP (rms or std -> no difference as mean is 0 in average reference)
    GFP_1 = std(topo_1,0,1);
    GFP_2 = std(topo_2,0,1);

    %norm topos
    topo_1 = topo_1./GFP_1; 
    topo_2 = topo_2./GFP_2; 

    %difference (rms or std: HERE IT MAKES A DIFFERENCE! mean of 
    %differences is not 0 anymore!!) Skrandies says that GMD is the std
    %between successive maps' difference!
    real_diff = std(topo_1 - topo_2,0,1);
    
    %% create a permutation distribution of differences
    %technically only a large bootstrap sample
    
    tempdata = cat(3, topos_1, topos_2);
    subs = size(tempdata,2);
    diff_distribution = nan(1,runs);
    
    for run = 1:runs
        subject_switch = randi([0 1], subs, 1);

        %switch condition means for subject
        for subject = find(subject_switch)        
            tempdata(:,subject,:) = flip(tempdata(:,subject,:),3);
        end

        %average, norm, subtract, std, store
        topo_1run = squeeze(mean(tempdata(:,:,1),2));
        topo_2run = squeeze(mean(tempdata(:,:,2),2));
        topo_1run = topo_1run./std(topo_1run,0,1);
        topo_2run = topo_2run./std(topo_2run,0,1);
        
        %basically simply the euclidean distance normed by n-1
        diff_distribution(1,run) = std(topo_1run - topo_2run,0,1);
    end
    
    %add real difference to distribution
    diff_distribution = [real_diff diff_distribution];
    
    %compute p-value
    try
        [f, x] = ecdf(diff_distribution);
        p_value = 1 - f(x == real_diff);
        cdf_available = 1;
    catch %if the statistics package is not installed
        p_value = sum(real_diff < diff_distribution)/size(diff_distribution,2);
        cdf_available = 10;
    end
    
    %logical output if topos are the same or different
    %NOTE: one-tailed t-test as we suspect the real difference to be larger
    %than most of the permutation distribution samples
    same = p_value > .05;

    %% plot result if requested

    if plot
        figure; 
        subplot(2,2,1); hold on; title('Topo 1 (normed)');
        topoplot(topo_1,chanlocs,...
            'electrodes','off','maplimits',maplim,'plotrad',plotrad,...
            'headrad',headrad,'intrad',intrad);
        axis([-0.6 0.6 -0.6 0.6]) 
        subplot(2,2,2); hold on; title('Topo 2 (normed)');
        topoplot(topo_2,chanlocs,...
            'electrodes','off','maplimits',maplim,'plotrad',plotrad,...
            'headrad',headrad,'intrad',intrad);
        axis([-0.6 0.6 -0.6 0.6])
        if cdf_available
            subplot(2,2,3); hold on;
            cdfplot(diff_distribution)
            line([real_diff real_diff],[0 1],'color',[1 0 0], 'linewidth',1)
        end
        subplot(2,2,4); hold on; title('Difference Distribution');
        hist(diff_distribution)
        ylim = get(gca,'ylim');        
        line([real_diff real_diff],ylim,'color',[1 0 0], 'linewidth',1)
    end
    
    return
