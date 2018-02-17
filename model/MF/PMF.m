function [ w1_P1, w1_M1 ] = PMF(train_vec, num_feat, num_p, num_m, maxepoch,...
    momentum, lambda, epsilon, batchsize, w1_P1, w1_M1, w1_P1_inc, w1_M1_inc )
if nargin < 6, momentum = 0.8; end
if nargin < 7, lambda = 0.1; end
if nargin < 8, epsilon = 0.5; end
if nargin < 9, batchsize = 1000; end
if nargin < 10, w1_P1 = 0.1*randn(num_p, num_feat); end
if nargin < 11, w1_M1 = 0.1*randn(num_m, num_feat); end
if nargin < 12, w1_P1_inc = zeros(num_p, num_feat); end
if nargin < 13, w1_M1_inc = zeros(num_m, num_feat); end

numdata = size(train_vec,1);
numbatches = ceil(numdata/batchsize);
for epoch = 1:maxepoch
    train_vec = train_vec(randperm(numdata),:);
    f_s = 0;
    %fprintf(1,'epoch %d\n',epoch);
    for batch = 1:numbatches
        aa_p   = double(train_vec((batch-1)*batchsize+1:min(batch*batchsize,numdata),1));
        aa_m   = double(train_vec((batch-1)*batchsize+1:min(batch*batchsize,numdata),2));
        rating = double(train_vec((batch-1)*batchsize+1:min(batch*batchsize,numdata),3));
        N = size(aa_p, 1);
        % rating = rating-mean_rating; % Default prediction is the mean rating.
        
        %%%%%%%%%%%%%% Compute Predictions %%%%%%%%%%%%%%%%%
        pred_out = sum(w1_M1(aa_m,:).*w1_P1(aa_p,:),2);
        %f = sum( (pred_out - rating).^2 + ...
        %    0.5*lambda*( sum( (w1_M1(aa_m,:).^2 + w1_P1(aa_p,:).^2),2)));
        
        %%%%%%%%%%%%%% Compute Gradients %%%%%%%%%%%%%%%%%%%
        IO = repmat(2*(pred_out - rating),1,num_feat);
        Ix_m=IO.*w1_P1(aa_p,:) + lambda*w1_M1(aa_m,:);
        Ix_p=IO.*w1_M1(aa_m,:) + lambda*w1_P1(aa_p,:);
        
        dw1_M1 = zeros(num_m,num_feat);
        dw1_P1 = zeros(num_p,num_feat);
        
        for ii=1:N
            dw1_M1(aa_m(ii),:) =  dw1_M1(aa_m(ii),:) +  Ix_m(ii,:);
            dw1_P1(aa_p(ii),:) =  dw1_P1(aa_p(ii),:) +  Ix_p(ii,:);
        end
        
        %%%% Update movie and user features %%%%%%%%%%%
        
        w1_M1_inc = momentum*w1_M1_inc + epsilon*dw1_M1/N;
        w1_M1 =  w1_M1 - w1_M1_inc;
        
        w1_P1_inc = momentum*w1_P1_inc + epsilon*dw1_P1/N;
        w1_P1 =  w1_P1 - w1_P1_inc;
        
        %%%%%%%%%%%%%% Compute Predictions after Paramete Updates %%%%%%%%%%%%%%%%%
        pred_out = sum(w1_M1(aa_m,:).*w1_P1(aa_p,:),2);
        f_s = f_s + sum( (pred_out - rating).^2 ); % + ...
           % 0.5*lambda*( sum( (w1_M1(aa_m,:).^2 + w1_P1(aa_p,:).^2),2)));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Compute predictions on the validation set %%%%%%%%%%%%%%%%%%%%%%
        % pred_out = sum(w1_M1(aa_m,:).*w1_P1(aa_p,:),2) + mean_rating;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    err_train = sqrt(f_s/numdata);
    fprintf(1, 'epoch %4i Training RMSE %6.4f\n', epoch, err_train);
end
end

