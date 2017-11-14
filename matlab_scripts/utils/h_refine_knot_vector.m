%% divide the knot vector V to num_division in each segments and do it num_times times
% V is assumed to be sorted ascending
function y = h_refine_knot_vector(V,num_division,num_times)

if(num_times > 0)
    Vk = unique(V);
    new_knots = zeros(1,(num_division-1)*(length(Vk)-1));
    k = 1;
    for i = 1:length(Vk)-1
        for j = 1:num_division-1
            new_knots(k) =  Vk(i) + (Vk(i+1)-Vk(i)) * j / (num_division);
            k = k + 1;
        end
    end
    y = sort([V new_knots]);
    y = h_refine_knot_vector(y,num_division,num_times-1);
else
    y = V;
end

