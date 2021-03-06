function M = predRustMT(p,s)
%M = predRustMT(p,s)
%
%Implements the MT model of Rust et al. (2006)

%s.dList is the list of possible stimulus directions
%p.pList is the list of preferred directions for each model V1 neuron

%% Step 1. 
%Calulate  the linear response to each possible direction based on tuning 
%curves for each neuron using the 'von Mises' function.  

%'d' will be a matrix with the rows are different neurons, and the columns 
%are the response of each neuron to each stimulus component. 
[thetam,pn] = meshgrid(s.dList,p.pList);
d = exp(p.b*cos( (thetam-pn)*pi/180));


%% Step 2. 
% normalize so the sum of each row is 1. Note: none of this has anything to
% do with the stimulus yet.
dprime = d./ repmat(sum(d,2),1,length(s.dList));

%%

%now, loop through each stimulus to get the model's response.
%Each row of s.c is a different stimulus, and contains a list of contrasts 
%for each of the directions.

for i=1:size(s.c,1)
    %Step 3. Calculate the linear response to the stimulus
    L = sum(dprime.*repmat(s.c(i,:),length(p.pList),1),2);
    %Step 4. Calculate the 'untuned' normalization
    P = L.^2 ./ (sum(L.^2) + p.s1^2);
    %Step 5. Calcualte the 'self' normalization
    V = P./(P+p.s2);
    %Step 6. Calculate the linearly weighted sum across V1 neurons
    Q = p.w*V;
    %Step 7. Calculate the nonlinear output
    M(i) = p.A.*exp(p.B*Q);
end

