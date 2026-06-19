%% Network Analysis of ADHD vs Control Participants
%This commented section is about the preprocessing of the raw data that was
%provided. The original data was a flattened data with a pairwise
%connection between the brain regions. The correlation used in getting the
%conenctions is Pearson. So basically I combined the flattened data and
%the clinical information data, and split by adhd_outcome into patient and
%control. So I uploaded the patient and control data instead of the
%original data.
%-----------------------------------------------------------------------

%reading the two datasets
%connectivity_data = readtable("Year2_Neuroscience/MATLAB/TRAIN_FUNCTIONAL_CONNECTOME_MATRICES_new_36P_Pearson.csv");
%adhd_labels = readtable("Year2_Neuroscience/MATLAB/adhd_labels.xlsx");
%Merging the data by participant_id
%merged_data = innerjoin(connectivity_data, adhd_labels, "Keys", "participant_id");
%dropping the Sex_F column
%merged_data.Sex_F = [];

%quick data quality checks
%connectivity_cols = merged_data{:, 2:end};
%nan_per_subject = sum(isnan(connectivity_cols), 2);
%fprintf('Subjects with missing values: %d out of %d\n', sum(nan_per_subject > 0), height(merged_data));
%fprintf('Max missing values in any subject: %d\n', max(nan_per_subject));
%nan_per_connection = sum(isnan(connectivity_cols), 1);
%fprintf('Connections with missing values: %d out of %d\n', sum(nan_per_connection > 0), size(connectivity_cols, 2));
%fprintf('Min value in data: %.4f\n', min(connectivity_cols(:), [], 'omitnan'));
%fprintf('Max value in data: %.4f\n', max(connectivity_cols(:), [], 'omitnan'))
%fprintf('Mean: %.4f, SD: %.4f, Median: %.4f\n', ...
%    mean(connectivity_cols(:)), std(connectivity_cols(:)), median(connectivity_cols(:)));
%fprintf('Perfect correlations: %.2f%%\n', 100*sum(connectivity_cols(:)==1)/numel(connectivity_cols));


%splitting the data by adhd_outcome into patient and control
%connectivity_data_patient = merged_data(merged_data.ADHD_Outcome==1, :);
%connectivity_data_control = merged_data(merged_data.ADHD_Outcome==0, :);

%converting connectivity to matrix data
%dropping the IDs and labels
%unwanted_columns = {'participant_id','ADHD_Outcome'};

%selecting only the connectivity features
%connectivity_values_ADHD = connectivity_data_patient(:, ...
%    ~ismember(connectivity_data_patient.Properties.VariableNames, unwanted_columns));

%connectivity_values_control = connectivity_data_control(:, ...
%    ~ismember(connectivity_data_control.Properties.VariableNames, unwanted_columns));

%converting tables to numeric matrices
%connectivity_matrix_ADHD = table2array(connectivity_values_ADHD);
%connectivity_matrix_control = table2array(connectivity_values_control);

%reconstructing the matrix
%function conn_matrix = reconstruct_matrix(flattened_vector)
%    no_of_regions = 200;
%    conn_matrix = zeros(no_of_regions, no_of_regions);
%   idx = 1;
%    for i = 1:no_of_regions
%        for j = (i+1):no_of_regions
%            conn_matrix(i,j) = flattened_vector(idx);
%            conn_matrix(j, i) = flattened_vector(idx);
%            idx = idx+1;
%        end
%    end
%end
%for patients
%no_of_patients = size(connectivity_matrix_ADHD, 1);
%fmri_ADHD = zeros(200, 200, no_of_patients);
%for i = 1:no_of_patients
%    fmri_ADHD(:,:,i) = reconstruct_matrix(connectivity_matrix_ADHD(i,:));
%end

%for controls
%no_of_controls = size(connectivity_matrix_control, 1);
%fmri_control = zeros(200, 200, no_of_controls);
%for i = 1:no_of_controls
%    fmri_control(:,:,i) = reconstruct_matrix(connectivity_matrix_control(i,:));
%end

%inspecting one subject per group visually
%imagesc(fmri_ADHD(:,:,6)); colorbar; title("ADHD Subject 1")

%figure; imagesc(fmri_control(:,:,6)); colorbar; title("Control Subject 1")

%sanity check for the no of subjects in each group
%size(fmri_control)
%size(fmri_ADHD)
%we have 
% ADHD: 831 subjects × 200 × 200 connectivity matrices
%Control: 382 subjects × 200 × 200 connectivity matrices
   
%saving the matrices
%save('ADHD_connectivity.mat','fmri_ADHD')
%save('Control_connectivity.mat','fmri_control')
%--------------------------------------------------------------------
%% Graph Analysis of ADHD vs Control Participants

load("ADHD_connectivity.mat", "fmri_ADHD");
load("Control_connectivity.mat", "fmri_control");

%setting the threshold
thresholds = [0.05 0.10 0.15 0.20 0.25 0.30];

%initializing arrays to store the metrics
num_ADHD = size(fmri_ADHD, 3);
num_control = size(fmri_control, 3);
no_of_permutations = 10000;
metrics = {'Pathlength', 'Clustering', 'GlobalEfficiency', 'LocalEfficiency', 'Smallworldness'};  
no_random_networks = 50;
no_rewire_iterations = 10;

for t = 1:length(thresholds)
    threshold = thresholds(t);
    fprintf('\n========================================\n');
    fprintf('=== Threshold = %.2f (%.0f%%) ===\n', threshold, threshold*100);
    fprintf('========================================\n');
    
    % Initialize arrays
    clustering_ADHD = zeros(num_ADHD, 1);
    path_length_ADHD = zeros(num_ADHD, 1);
    globaleff_ADHD = zeros(num_ADHD, 1);
    localeff_ADHD = zeros(num_ADHD, 1);
    smallworld_ADHD = zeros(num_ADHD, 1);

    clustering_control = zeros(num_control, 1);
    path_length_control = zeros(num_control, 1);
    globaleff_control = zeros(num_control, 1);
    localeff_control = zeros(num_control, 1);
    smallworld_control = zeros(num_control, 1);

    % GENERATE RANDOM NETWORKS ONCE FOR THIS THRESHOLD
    fprintf('\nGenerating %d random networks for comparison...\n', no_random_networks);
    %the idea of using a subject as template to generate random network for small worldness calculation was by claude
    %given that my previous attempt at generating random network per
    %subject and threshold took hours to compute. This method reduced the
    %computational time and made processing faster.


    % Use first ADHD subject as template
    W_template = fmri_ADHD(:,:,1);
    W_template(W_template<0) = 0;
    W_template = threshold_proportional(W_template, threshold);
    
    % Generate random networks
    clustering_random = zeros(no_random_networks, 1);
    pathlength_random = zeros(no_random_networks, 1);
    
    for r = 1:no_random_networks
        W_rand = randmio_und_connected(W_template, no_rewire_iterations);
        L_rand = weight_conversion(W_rand, 'lengths');
        D_rand = distance_wei(L_rand);
        
        clustering_random(r) = mean(clustering_coef_wu(W_rand));
        pathlength_random(r) = charpath(D_rand, 0, 0);
        
        if mod(r, 10) == 0
            fprintf('  Generated %d/%d random networks\n', r, no_random_networks);
        end
    end
    
    mean_clustering_random = mean(clustering_random);
    mean_pathlength_random = mean(pathlength_random);
    
    fprintf('Random network metrics: C_rand = %.4f, L_rand = %.4f\n', ...
        mean_clustering_random, mean_pathlength_random);
    
    % COMPUTE METRICS FOR ADHD GROUP
    fprintf('\nProcessing ADHD subjects...\n');
    
    for subj = 1:num_ADHD
        W = fmri_ADHD(:,:,subj);
        W(W<0) = 0;
        W = threshold_proportional(W, threshold);
        L = weight_conversion(W, 'lengths');
    
        clustering_per_node = clustering_coef_wu(W);
        D = distance_wei(L);
        path_length = charpath(D, 0, 0);
        globalEff = efficiency_wei(W);
        localEff = mean(efficiency_wei(W,1));
        
        clustering_ADHD(subj) = mean(clustering_per_node);
        path_length_ADHD(subj) = path_length;
        globaleff_ADHD(subj) = globalEff;
        localeff_ADHD(subj) = localEff;
        
        % Compute small-worldness using pre-computed random metrics
        gamma = clustering_ADHD(subj) / mean_clustering_random;
        lambda = path_length_ADHD(subj) / mean_pathlength_random;
        smallworld_ADHD(subj) = gamma / lambda;
        
        if mod(subj, 100) == 0
            fprintf('  Completed %d/%d ADHD subjects\n', subj, num_ADHD);
        end
    end
    
    % COMPUTE METRICS FOR CONTROL GROUP
    fprintf('\nProcessing Control subjects...\n');
    
    for subj = 1:num_control
        W = fmri_control(:,:,subj);
        W(W<0) = 0;
        W = threshold_proportional(W, threshold);
        L = weight_conversion(W, 'lengths');

        clustering_per_node = clustering_coef_wu(W);
        D = distance_wei(L);
        path_length = charpath(D, 0, 0);
        globalEff = efficiency_wei(W);
        localEff = mean(efficiency_wei(W,1));
        
        clustering_control(subj) = mean(clustering_per_node);
        path_length_control(subj) = path_length;
        globaleff_control(subj) = globalEff;
        localeff_control(subj) = localEff;
        
        % Compute small-worldness using pre-computed random metrics
        gamma = clustering_control(subj) / mean_clustering_random;
        lambda = path_length_control(subj) / mean_pathlength_random;
        smallworld_control(subj) = gamma / lambda;
        
        if mod(subj, 50) == 0
            fprintf('  Completed %d/%d Control subjects\n', subj, num_control);
        end
    end

    % REGRESS OUT MEAN FC
    fprintf('\nRegressing out mean functional connectivity...\n');
    
    % Compute mean FC for each subject
    mean_FC_ADHD = zeros(num_ADHD, 1);
    mean_FC_control = zeros(num_control, 1);
    
    for subj = 1:num_ADHD
        W = fmri_ADHD(:,:,subj);
        W(W<0) = 0;
        upper_tri_idx = triu(true(size(W)), 1);
        mean_FC_ADHD(subj) = mean(W(upper_tri_idx));
    end
    
    for subj = 1:num_control
        W = fmri_control(:,:,subj);
        W(W<0) = 0;
        upper_tri_idx = triu(true(size(W)), 1);
        mean_FC_control(subj) = mean(W(upper_tri_idx));
    end
    
    % Store original metrics for comparison
    path_length_ADHD_original = path_length_ADHD;
    clustering_ADHD_original = clustering_ADHD;
    globaleff_ADHD_original = globaleff_ADHD;
    localeff_ADHD_original = localeff_ADHD;
    smallworld_ADHD_original = smallworld_ADHD;

    path_length_control_original = path_length_control;
    clustering_control_original = clustering_control;
    globaleff_control_original = globaleff_control;
    localeff_control_original = localeff_control;
    smallworld_control_original = smallworld_control;
    
    % Regress out mean FC from ADHD metrics
    %The part of regressing out the mean from each group was generated by
    %claude
    
    X_ADHD = [ones(num_ADHD, 1), mean_FC_ADHD];
    
    beta = X_ADHD \ clustering_ADHD;
    residuals = clustering_ADHD - X_ADHD * beta;
    clustering_ADHD = residuals + mean(clustering_ADHD);
    
    beta = X_ADHD \ path_length_ADHD;
    residuals = path_length_ADHD - X_ADHD * beta;
    path_length_ADHD = residuals + mean(path_length_ADHD);
    
    beta = X_ADHD \ globaleff_ADHD;
    residuals = globaleff_ADHD - X_ADHD * beta;
    globaleff_ADHD = residuals + mean(globaleff_ADHD);
    
    beta = X_ADHD \ localeff_ADHD;
    residuals = localeff_ADHD - X_ADHD * beta;
    localeff_ADHD = residuals + mean(localeff_ADHD);
    
    beta = X_ADHD \ smallworld_ADHD;
    residuals = smallworld_ADHD - X_ADHD * beta;
    smallworld_ADHD = residuals + mean(smallworld_ADHD);
    
    % Regress out mean FC from Control metrics
    X_control = [ones(num_control, 1), mean_FC_control];
    
    beta = X_control \ clustering_control;
    residuals = clustering_control - X_control * beta;
    clustering_control = residuals + mean(clustering_control);
    
    beta = X_control \ path_length_control;
    residuals = path_length_control - X_control * beta;
    path_length_control = residuals + mean(path_length_control);
    
    beta = X_control \ globaleff_control;
    residuals = globaleff_control - X_control * beta;
    globaleff_control = residuals + mean(globaleff_control);
    
    beta = X_control \ localeff_control;
    residuals = localeff_control - X_control * beta;
    localeff_control = residuals + mean(localeff_control);

    beta = X_control \ smallworld_control;
    residuals = smallworld_control - X_control * beta;
    smallworld_control = residuals + mean(smallworld_control);
    
    % PERMUTATION TESTS WITH REGRESSED MEAN FC
    fprintf('\n--- Results WITH regressed mean FC ---\n');
    ADHD_data = {path_length_ADHD, clustering_ADHD, globaleff_ADHD, localeff_ADHD, smallworld_ADHD};
    control_data = {path_length_control, clustering_control, globaleff_control, localeff_control, smallworld_control};
    run_permutation_test(ADHD_data, control_data, no_of_permutations, metrics);

    % PERMUTATION TESTS WITHOUT REGRESSION
    fprintf('\n--- Results WITHOUT regressed mean FC ---\n');
    ADHD_data_orig = {path_length_ADHD_original, clustering_ADHD_original, globaleff_ADHD_original, localeff_ADHD_original, smallworld_ADHD_original};
    control_data_orig = {path_length_control_original, clustering_control_original, globaleff_control_original, localeff_control_original, smallworld_control_original};
    run_permutation_test(ADHD_data_orig, control_data_orig, no_of_permutations, metrics);
    
    % SMALL-WORLDNESS SUMMARY
    fprintf('\n--- Small-Worldness Summary ---\n');
    
    % Count subjects with small-world properties (sigma > 1)
    smallworld_threshold = 1.0;
    num_sw_ADHD = sum(smallworld_ADHD_original > smallworld_threshold);
    num_sw_control = sum(smallworld_control_original > smallworld_threshold);
    
    pct_sw_ADHD = 100 * num_sw_ADHD / num_ADHD;
    pct_sw_control = 100 * num_sw_control / num_control;
    
    fprintf('Networks with small-world properties (σ > %.1f):\n', smallworld_threshold);
    fprintf('  ADHD: %.1f%% (%d/%d subjects)\n', pct_sw_ADHD, num_sw_ADHD, num_ADHD);
    fprintf('  Control: %.1f%% (%d/%d subjects)\n', pct_sw_control, num_sw_control, num_control);
    
    fprintf('\nMean small-worldness (σ):\n');
    fprintf('  ADHD: %.3f (SD = %.3f)\n', mean(smallworld_ADHD_original), std(smallworld_ADHD_original));
    fprintf('  Control: %.3f (SD = %.3f)\n', mean(smallworld_control_original), std(smallworld_control_original));
    
    fprintf('\nMedian small-worldness (σ):\n');
    fprintf('  ADHD: %.3f\n', median(smallworld_ADHD_original));
    fprintf('  Control: %.3f\n', median(smallworld_control_original));
end

fprintf('\n========================================\n');
fprintf('Analysis Complete!\n');
fprintf('========================================\n');

% Permutation test Function
function [observed_difference, p_values] = run_permutation_test(ADHD_data, control_data, no_of_permutations, metrics)
    n_ADHD = length(ADHD_data{1});
    n_control = length(control_data{1});
   
    for m = 1:length(metrics)
        observed_difference = mean(ADHD_data{m}) - mean(control_data{m});
        combined_data = [ADHD_data{m}; control_data{m}];
        permutation_difference = zeros(no_of_permutations, 1);
        
        for p = 1:no_of_permutations
            permuted_labels = combined_data(randperm(length(combined_data)));
            permutation_ADHD = permuted_labels(1:n_ADHD);
            permutation_control = permuted_labels(n_ADHD+1:end);
            permutation_difference(p) = mean(permutation_ADHD) - mean(permutation_control);
        end
        
        % Two-tailed p-value
        p_values = mean(abs(permutation_difference) >= abs(observed_difference));
        
        fprintf('%s: Observed difference = %.4f, p-value = %.4f\n', ...
            metrics{m}, observed_difference, p_values);
    end
end