% shShowMtPopulationResponse(populationResponse, vx, vy)
%
% Show the response of a population of MT neurons that was computed using
% shMtPopulationResponse.
%
% Required arguments:
% populationResponse            the output from shMtPopulationResponse
% vx                            the second output from shMtPopulationResponse
% vx                            the third output from shMtPopulationResponse

function shShowMtPopulationResponse(populationResponse, vx, vy);

if size(vx, 1) > 1 & size(vx, 2)
    vx = vx(1,:);
    vy = vy(:,1);
end
image(vx, vy, gray2rgbsc(populationResponse));

xlabel('V_x');
ylabel('V_y');