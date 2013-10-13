function [digit_mats] = make_digit_mats(filename)

fid = fopen(filename);
counter = zeros(10);

while 1

	tline = fgetl(fid);
	if ~ischar(tline)
		break;
	end

	place_holder = sscanf(tline, '%f');
	dig_num = place_holder(1, 1);

	if dig_num == 0
		dig_num = 10;
	end

	counter(dig_num) = counter(dig_num) + 1;
	digit_mats(:, counter(dig_num), dig_num) = place_holder(2:257);
end

fclose(fid);

digit_mats = 0.5.*(digit_mats .+ 1);