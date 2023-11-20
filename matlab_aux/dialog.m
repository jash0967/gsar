function [x] = dialog()
%DIALOG Summary of this function goes here
%   Detailed explanation goes here
prompt = "Choose the variant: \n 1. Pure complementarity \n 2. Pure conformity \n 3. Generalized SAR \n >> ";
x=input(prompt);

disp("List of the networks: \n 1. Money \n 2. Advice \n 3. Kerosine \n 4. Medical \n 5. Non-relatives \n 6. Relatives \n 7. Temple \n 8. Visit")
prompt = "Choose the networks. All will be used if empty. \n >> ";
x=input(prompt);
if x==empty
    x=1:8;
else
    
end


end
