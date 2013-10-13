function [stats, success_rate] = run_algorithm(algo_type, digit_mats)

rnk = 10;

if ~exist('digit_mats')
	digit_mats = make_digit_mats('data/zip.train');
end

if strcmp(algo_type, 'svd')
	for i = 1:10
		[u, s, v] = svd(digit_mats(:, :, i));
		M(:, :, i) = u(:, 1:rnk);
	end
elseif strcmp(algo_type, 'nmf')
	for i = 1:10
		[w, h] = nmf_mu(digit_mats(:, :, i), rnk, 100);
		M(:, :, i) = w;
	end
end

fid = fopen('data/zip.test');

stats = zeros(10, 5); stats(:, 1) = 1:10; stats(10, 1) = 0;

while 1
	
	tline = fgetl(fid);
	if ~ischar(tline)
		break;
	end

	td = sscanf(tline, '%f');

	digit_num = td(1, 1);
	if digit_num == 0
		digit_num = 10;
	end

	test_dig = td(2:257);
	test_dig = 0.5.*(test_dig .+ 1);

	stats(digit_num, 2) = stats(digit_num, 2) + 1;

	resids = zeros(10, 1);

	if strcmp(algo_type, 'svd')
		for i = 1:10
			x = M(:, :, i)'*test_dig;
			resids(i) = norm(M(:, :, i)*x - test_dig);
		end
	elseif strcmp(algo_type, 'nmf')
		for i = 1:10
			x = M(:, :, i)\test_dig;
			resids(i) = norm(M(:, :, i)*x - test_dig);
		end
	end

	[smallest_resid, digit_result] = min(resids);

	if digit_num == digit_result
		stats(digit_num, 3) = stats(digit_num, 3) + 1;
	else
		stats(digit_num, 4) = stats(digit_num, 4) + 1;
	end
end

fclose(fid);

stats(:, 5) = stats(:, 3)./stats(:, 2);
success_rate = sum(stats(:, 3))/sum(stats(:, 2));
