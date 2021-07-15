function [vq] = Interp2_TSE(x, y, v, xq, yq, N, relNum, method)
%INTERP2_TSE interpolates 2d data using the Taylor Series Expansion based method
%
%   Notes:
%         1) if the distance between the query point and the given point is less
%            than 0.1, the program treat them as the same point. The threshold is 
%            related to the positioning accuracy. 0.1 is more than enough for most applications.
%         2) in current version, we limit the value of 'N' less than 9.
%         3) To improve the stability of the numerical solution, we assemble a over-determinated
%            matrix equation, where, the number of eqs is 120% of the unknown in default.
%         4) when the terrain varies greatly. we recommend setting method to 'nearest'
%         5) x/y/v/xq/yq can be a vector.
%
%   Input:
%         x                   required, matrix, double      coordinates of original data in the north
%         y                   required, matrix, double      coordinates of original data in the east
%         v                   required, matrix, double      value of original data
%         xq                  required, matrix, double      coordinates of interpolated data in the north
%         yq                  required, matrix, double      coordinates of interpolated data in the east
%         N                   required, scalar, double      order of the Taylor series
%         relNum              optional, scalar, double      num of eqs relative to the num of ders
%         method              optional, vector, char        ['uniform' 'nearest'], method to select neighbors of the given point
%   Output:
%         vq                  matrix, double                value of interpolated data
%
%   Reference:
%         2021, Chen and Yang, Potential field data reconstruction using Taylor series
%               expansion.
%
%   Update records:
%         author        date             features
%         Tao Chen      03-Nov-2020      standard implementation, improve stability, parallelization
%         Tao Chen      11-Mar-2021      add detailed description
%
%   See also Interp1_TSE, Interp3_TSE.
%
%   Copyright Tao Chen.

%% Parse and initialize the optional parameters
narginchk(5, 8);
if 5 == nargin
    N = 5;
    relNum = 1.2;
    method = 'uniform';
elseif 6 == nargin
    relNum = 1.2;
    method = 'uniform';
elseif 7 == nargin
    method = 'uniform';
end

if N > 9
    error('Interp2_TSE: FunctionInput: Invalid N (N <= 9)');
end

N = N + 1;   % we include the remainder of the Taylor series when assembling the matrix equation

% convert input variables to column-vectors
[nr, nc] = size(xq);
x = x(:);
y = y(:);
v = v(:);
xq = xq(:);
yq = yq(:);
vq = zeros(size(xq));   % initialization

% coefficients of the matrix equation, details refer function 'CoeffsOfTSE2D'. For saving time, we pre-save coefficients with order less than 9
tseCoeffs = [1,0,1;0,1,1;2,0,0.500000000000000;1,1,1;0,2,0.500000000000000;3,0,0.166666666666667;2,1,0.500000000000000;1,2,0.500000000000000;0,3,0.166666666666667;4,0,0.0416666666666667;3,1,0.166666666666667;2,2,0.250000000000000;1,3,0.166666666666667;0,4,0.0416666666666667;5,0,0.00833333333333333;4,1,0.0416666666666667;3,2,0.0833333333333333;2,3,0.0833333333333333;1,4,0.0416666666666667;0,5,0.00833333333333333;6,0,0.00138888888888889;5,1,0.00833333333333333;4,2,0.0208333333333333;3,3,0.0277777777777778;2,4,0.0208333333333333;1,5,0.00833333333333333;0,6,0.00138888888888889;7,0,0.000198412698412698;6,1,0.00138888888888889;5,2,0.00416666666666667;4,3,0.00694444444444444;3,4,0.00694444444444444;2,5,0.00416666666666667;1,6,0.00138888888888889;0,7,0.000198412698412698;8,0,2.48015873015873e-05;7,1,0.000198412698412698;6,2,0.000694444444444445;5,3,0.00138888888888889;4,4,0.00173611111111111;3,5,0.00138888888888889;2,6,0.000694444444444445;1,7,0.000198412698412698;0,8,2.48015873015873e-05;9,0,2.75573192239859e-06;8,1,2.48015873015873e-05;7,2,9.92063492063492e-05;6,3,0.000231481481481482;5,4,0.000347222222222222;4,5,0.000347222222222222;3,6,0.000231481481481482;2,7,9.92063492063492e-05;1,8,2.48015873015873e-05;0,9,2.75573192239859e-06;10,0,2.75573192239859e-07;9,1,2.75573192239859e-06;8,2,1.24007936507937e-05;7,3,3.30687830687831e-05;6,4,5.78703703703704e-05;5,5,6.94444444444444e-05;4,6,5.78703703703704e-05;3,7,3.30687830687831e-05;2,8,1.24007936507937e-05;1,9,2.75573192239859e-06;0,10,2.75573192239859e-07];
% tseCoeffs = CoeffsOfTSE2D(N);

numDers = sum((1 : N) + 1);   % number of the derivatives
numEqs = floor(numDers * relNum);
numDersForRec = sum((1 : N - 1) + 1);   % the remainder of the Taylor series is not used for the reconstruction

% interpolation
parfor i = 1 : length(xq)

    [valMin, indMin] = min(sqrt((xq(i) - x).^2 + (yq(i) - y).^2));   % the nearest neighbor of the query point is set as given point
    
    if abs(valMin) <= 1e-1   % the query point is actually the given point, 
        vq(i) = v(indMin);
    else

        % 1-1) select numEqs neighbors of the given point to assemble the coefficients of matrix equation
        %    make the selected neighbors of the given point evenly distributed in all directions
        if strcmpi(method, 'uniform')
            neighborhoodsIndex = setdiff(1 : length(x), indMin);   % all available neighbors of the given point
            [theta, r] = cart2pol(y(neighborhoodsIndex) - y(indMin), x(neighborhoodsIndex) - x(indMin));   % locations of all neighborhoods in polar coordinates
            theta = rad2deg(theta);
            
            thetaRangeOfSections = floor(min(theta)) : (ceil(max(theta)) - floor(min(theta)))/numEqs : ceil(max(theta));   % divide the neighborhoods into 'numEqs' sections
            candNeighborhoodsIndex = zeros(numEqs, numEqs);
            
            for jj = 1 : numEqs   % loop for sections
                
                neighborhoodsIndOfSections = find(theta >= thetaRangeOfSections(jj) & theta < thetaRangeOfSections(jj + 1));   % neighborhoods for current section
                
                % 'numEqs' candidate neighborhoods of current section
                [~, candNeighborhoodsIndOfSections] = mink(r(neighborhoodsIndOfSections), numEqs);   % the candidates near to the given point
                %             [~, candNeighborhoodsIndOfSections] = mink(abs(r(neighborhoodsIndOfSections) - valMin), numEqs);   % the distance between the given point and candidates are close to 'valMin'
                %             [~, candNeighborhoodsIndOfSections] = mink(abs(round(r(neighborhoodsIndOfSections)/valMin) - r(neighborhoodsIndOfSections)/valMin), numEqs);   % candidates shall the same space resolution
                
                if ~isempty(candNeighborhoodsIndOfSections)
                    candNeighborhoodsIndex(1 : length(candNeighborhoodsIndOfSections), jj) = ...
                        neighborhoodsIndex(neighborhoodsIndOfSections(candNeighborhoodsIndOfSections));   % all candidata neighborhoods for all sections
                end
                
            end
            indCand = candNeighborhoodsIndex';
            indCand = nonzeros(indCand);
            ind = indCand(1 : numEqs);
        end
        
        % 1-2) nearest princple
        if strcmpi(method, 'nearest')
            [~, r] = cart2pol(y - y(indMin), x - x(indMin));
            [~, ind] = mink(r, numEqs + 1);
            ind = setdiff(ind, indMin);
        end
        
        % plot function to visualize the distributation of the selected neighborhoods.
%         if i == 200 || i == 4446   % set the number what you want
%             figure;
%             axis off;
%             scatter(y, x, 50, 'k', 'filled'); hold on;
%             scatter(yq(i), xq(i), 80, 'g', 'filled'); 
%             scatter(y(indMin), x(indMin), 80, [1 0 1], 'filled');            
%             scatter(y(ind), x(ind), 80, 'b', 'filled');
%             
%             for jj = 1 : numEqs 
%                 [east_sum, north_sum] = pol2cart(deg2rad(thetaRangeOfSections(jj:jj+1)), 50);
%                 for kk = 1 : length(north_sum)
%                     plot([y(indMin) y(indMin)+east_sum(kk)], [x(indMin) x(indMin)+north_sum(kk)], '-r', 'LineWidth', 2.5);
%                 end
%             end
%         end
        
        detX = x(ind) - x(indMin);
        detY = y(ind) - y(indMin);
        rhsOfEqs = v(ind) - v(indMin);

        % 2) assemble the coefficients of the matrix equation
        coeffMtx = detX(:, ones(1, numDers)) .^ (tseCoeffs(1 : numDers, ones(1, numEqs))') .* ...
                   detY(:, ones(1, numDers)) .^ (tseCoeffs(1 : numDers, 2 * ones(1, numEqs))') .* ...
                   (tseCoeffs(1 : numDers, 3 * ones(1, numEqs))');
               
               
        % 3) get derivatives of the given point polyfit
        [fsOfEqs, ~] = linsolve(coeffMtx, rhsOfEqs);
%         [fsOfEqs, ~] = lsqr(coeffMtx, rhsOfEqs);
        
        % 4) reconstruction of value of query point
        detX = xq(i) - x(indMin);
        detY = yq(i) - y(indMin);

        vq(i) = (detX .^ tseCoeffs(1 : numDersForRec, 1)') .* ...
                (detY .^ tseCoeffs(1 : numDersForRec, 2)') .* ...
                tseCoeffs(1 : numDersForRec, 3)' * ...
                fsOfEqs(1 : numDersForRec) + ...
                v(indMin);
        
    end

end

vq = reshape(vq, nr, nc);

end

function tseCoeffs = CoeffsOfTSE2D(N)   %#ok<*DEFNU>
%COEFFSOFTSE2D generates the coefficients of the taylor series of order N (detX^b * detY^c * d)
%
%   Input:
%         N               required, scalar, double      order of the TSE
%   Output:
%         tseCoeffs       required, matrix, double      the sequence is [b, c, d]

tseCoeffs = zeros(N, 3);
for iN = 1 : N
    for i = iN : -1 : 0   % x    
        for j = iN - i : -1 : 0   % y
            if i + j == iN
                temp = [ones(1, i), 2 * ones(1, j)];
                tseCoeffs(iN, :) = [i, j, size(unique(perms(temp), 'row'), 1) / factorial(iN)];
            end        
        end
    end
end

end