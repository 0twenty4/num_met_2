#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

enum METHOD {
	JACOBI_RELAX_PARAM,
	GAUSS_ZEIDEL_RELAX_PARAM,
};

void input_diag(double**& sys_mat, double*& sys_right_vec, double*& sys_sol_begin_aprox,
	int& n, int& m, double& epsilon, int& max_iter,
	std::string sys_params_file, std::string sys_mat_file, std::string sys_right_vec_file,
	std::string sys_sol_begin_aprox_file) {

	std::ifstream in(sys_params_file);

	in >> n >> m >> epsilon >> max_iter;

	sys_mat = new double* [n];

	in.close();

	in.open(sys_mat_file);

	while (!in.eof()) {

		int first_col = 5;
		int shift = 0;
		int skip_elem;
		int source_mat_row_elem;
		for (int i = 0; i < n; i++) {

			sys_mat[i] = new double[9] {};

			source_mat_row_elem = 0;

			if (i == m + 2)
				shift = 0;

			if (i > 1 && i < m + 2 || i>m + 4)
				shift++;

			if (shift == 0)
				first_col--;

			for (int skipped_elems = 0; skipped_elems < shift; skipped_elems++) {
				in >> skip_elem;
				source_mat_row_elem++;
			}

			for (int j = first_col; source_mat_row_elem < n; j++) {
				int temp; in >> temp;
				sys_mat[i][j] = temp;

				source_mat_row_elem++;

				if (j == 2 or j == 5 or j == 8) {
					for (int skipped_elems = 0; skipped_elems < m and source_mat_row_elem < n; skipped_elems++) {
						in >> skip_elem;
						source_mat_row_elem++;
					}
				}
			}
		}
	}



	in.close();
	in.open(sys_right_vec_file);

	sys_right_vec = new double[n];

	for (int i = 0; i < n; i++) {
		in >> sys_right_vec[i];
	}

	in.close();
	in.open(sys_sol_begin_aprox_file);

	sys_sol_begin_aprox = new double[n];

	for (int i = 0; i < n; i++) {
		in >> sys_sol_begin_aprox[i];
	}

	in.close();
}

void input_block(double****& sys_mat, int& n, double**& sys_right_vec, double**& sys_sol_begin_aprox,
	double& epsilon, int& max_iter, int block_dim,
	std::string sys_params_file, std::string sys_mat_file, std::string sys_right_vec_file,
	std::string sys_sol_begin_aprox_file) {

	std::ifstream in(sys_params_file);
	in >> n >> epsilon >> max_iter;
	in.close();

	int block_size = block_dim * block_dim;
	int block_mat_dim = n / block_dim;

	in.open(sys_mat_file);

	sys_mat = new double*** [block_mat_dim];

	int block;
	while (!in.eof()) {
		for (int block_i = 0; block_i < block_mat_dim; block_i++) {
			sys_mat[block_i] = new double** [block_mat_dim];
			for (int i = 0; i < block_dim; i++) {
				for (int block_j = 0; block_j < block_mat_dim; block_j++) {

					if (i == 0)
						sys_mat[block_i][block_j] = new double* [block_dim];

					sys_mat[block_i][block_j][i] = new double[block_dim];

					for (int j = 0; j < block_dim; j++) {
						int temp; in >> temp;
						sys_mat[block_i][block_j][i][j] = temp;
					}

				}
			}
		}
	}

	in.close();
	in.open(sys_right_vec_file);

	sys_right_vec = new double* [block_mat_dim];

	for (int block_i = 0; block_i < block_mat_dim; block_i++) {
		sys_right_vec[block_i] = new double[block_dim];

		for (int i = 0; i < block_dim; i++) {
			in >> sys_right_vec[block_i][i];
		}
	}

	in.close();
	in.open(sys_sol_begin_aprox_file);

	sys_sol_begin_aprox = new double* [block_mat_dim];

	for (int block_i = 0; block_i < block_mat_dim; block_i++) {
		sys_sol_begin_aprox[block_i] = new double[block_dim];

		for (int i = 0; i < block_dim; i++) {
			in >> sys_sol_begin_aprox[block_i][i];
		}
	}

	in.close();

}

void output_vec(double* vec, int n, std::string output_file) {
	std::ofstream out(output_file);

	out << std::fixed << std::setprecision(15);

	for (int i = 0; i < n; i++)
		out << vec[i] << std::endl;
}

void output_block_vec(double** block_vec, int n, int block_dim, std::string output_file) {
	std::ofstream out(output_file);

	out << std::fixed << std::setprecision(15);

	int block_vec_dim = n / block_dim;

	for (int block_i = 0; block_i < block_vec_dim; block_i++) {
		for (int i = 0; i < block_dim; i++) {
			out << block_vec[block_i][i] << std::endl;
		}
	}
}

double* mult_mat_diag_vec(double** mat, double* vec, int n, int m) {

	double* res = new double[n] {};

	for (int i = 0, shift = 0; i < n; i++) {
		// достигнут нижний блок ненулевых диагоналей
		if (i == m + 2)
			shift = 0;

		// строки, включающие значения нулевых диагоналей между стредним и нижним блоком или ниже нижнего блока
		if (i > 1 && i < m + 2 || i>m + 4)
			shift++;

		for (int j = 0, vec_i = 0; j < n; j++) {
			if (mat[i][j] != 0) {

				res[i] += mat[i][j] * vec[vec_i + shift];

				vec_i++;

				//переход от одного блока к другому
				if ((j + 1) % 3 == 0)
					vec_i += m;

			}
		}
	}

	return res;

}

double* vec_substr(double* vec1, double* vec2, int n) {
	double* res = new double[n];

	for (int i = 0; i < n; i++) {
		res[i] = vec1[i] - vec2[i];
	}

	return res;
}

double* vec_sum(double* vec1_elems, double* vec2_elems, int n) {
	double* res = new double[n];

	for (int i = 0; i < n; i++) {
		res[i] = vec1_elems[i] + vec2_elems[i];
	}

	return res;
}

double vec_norm(double* vec, int n) {

	double sum = 0.;

	for (int i = 0; i < n; i++) {
		sum += vec[i] * vec[i];
	}

	return sqrt(sum);
}

double scal_mult(double* vec1, double* vec2, int n) {
	double res = 0.;

	for (int i = 0; i < n; i++) {
		res += vec1[i] * vec2[i];
	}

	return res;
}

double scal_mult_diag(double* vec1_diag, int i, double* vec2, int n, int m) {
	int k = 9;
	int shift = 0;
	double res = 0;

	if (i > 1 and i < 2 + m)
		shift = i - 1;
	else if (i > m + 4)
		shift = i - (m + 4);

	for (int j = 0, vec_j = 0; j < k; j++) {
		if (vec1_diag[j] != 0) {

			res += vec1_diag[j] * vec2[vec_j + shift];

			vec_j++;

			//переход от одного блока к другому
			if ((j + 1) % 3 == 0)
				vec_j += m;

		}
	}

	return res;
}

double* block_mult_mat_row_vec(double*** mat_row, double** vec, int n, int block_dim) {
	int block_mat_dim = n / block_dim;

	double* res = new double[block_dim] {};

	for (int block_j = 0; block_j < block_mat_dim; block_j++) {
		for (int i = 0; i < block_dim; i++) {
			res[i] += scal_mult(mat_row[block_j][i], vec[block_j], block_dim);
		}
	}

	return res;
}

double* mult_num_vec(double num, double* vec, int n) {
	double* mult = new double[n];
	for (int i = 0; i < n; i++) {
		mult[i] = num * vec[i];
	}

	return mult;
}


double relative_discrepancy(double** sys_mat, int n, int m, double* sys_right_vec, double* sys_sol) {
	return vec_norm(vec_substr(sys_right_vec, mult_mat_diag_vec(sys_mat, sys_sol, n, m), n), n) / vec_norm(sys_right_vec, n);
}

double* mult_mat_vec(double** mat, double* vec, int n) {
	double* mult = new double[n] {};

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			mult[i] += mat[i][j] * vec[j];
		}
	}

	return mult;
}

void gen_conditional(int n, double* vec, int k,
	double* not_diag_elems, int not_diag_elems_size,
	int rand_seed, std::string sys_mat_file, std::string sys_right_vec_file) {

	std::srand(rand_seed);

	double** sys_mat = new double* [n];

	for (int i = 0; i < n; i++) {
		sys_mat[i] = new double[n] {};
		for (int j = 0; j < n; j++) {

			if (i != j) {
				sys_mat[i][j] = not_diag_elems[std::rand() % not_diag_elems_size];
				sys_mat[i][i] -= sys_mat[i][j];
			}
		}

	}

	sys_mat[0][0] += 1;

	double* sys_right_vec = mult_mat_vec(sys_mat, vec, n);

	std::ofstream out(sys_mat_file);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			out << sys_mat[i][j] << " ";
		}
		out << std::endl;
	}

	out.close();
	out.open(sys_right_vec_file);

	for (int i = 0; i < n; i++) {
		out << sys_right_vec[i] << std::endl;
	}

}

double* mult_mat_LDU_vec(double** mat, double* vec, int n) {
	double* inter_res = new double[n] {};

	for (int i = 0; i < n; i++) {
		inter_res[i] += vec[i];
		for (int j = i + 1; j < n; j++) {
			inter_res[i] += mat[i][j] * vec[j];

		}
	}

	for (int i = 0; i < n; i++)
		inter_res[i] *= mat[i][i];

	double* res = new double[n] {};
	for (int i = 0; i < n; i++) {
		res[i] += inter_res[i];
		for (int j = 0; j < i; j++) {
			res[i] += mat[i][j] * inter_res[j];

		}
	}
	return res;
}

double block_relative_discrepancy(double**** sys_mat, int n, int block_dim, double** sys_right_vec, double** sys_sol) {
	int block_mat_dim = n / block_dim;
	double** mult_sys_mat_sys_sol = new double* [block_mat_dim];

	for (int block_i = 0; block_i < block_mat_dim; block_i++) {
		mult_sys_mat_sys_sol[block_i] = new double[block_dim] {};
		for (int block_j = 0; block_j < block_mat_dim; block_j++) {
			if (block_i == block_j)
				mult_sys_mat_sys_sol[block_i] = vec_sum(mult_sys_mat_sys_sol[block_i], mult_mat_LDU_vec(sys_mat[block_i][block_i], sys_sol[block_i], block_dim), block_dim);
			else {
				mult_sys_mat_sys_sol[block_i] = vec_sum(mult_sys_mat_sys_sol[block_i], mult_mat_vec(sys_mat[block_i][block_j], sys_sol[block_j], block_dim), block_dim);
			}
		}
	}

	double* substr;
	double sum1 = 0;
	double sum2 = 0;
	for (int block_i = 0; block_i < block_mat_dim; block_i++) {
		substr = vec_substr(sys_right_vec[block_i], mult_sys_mat_sys_sol[block_i], block_dim);
		sum1 += scal_mult(substr, substr, block_dim);
		sum2 += scal_mult(sys_right_vec[block_i], sys_right_vec[block_i], block_dim);
	}

	return sqrt(sum1) / sqrt(sum2);


}

void iter_step(METHOD method, double** sys_mat, int n, int m, double* sys_right_vec, double* sys_sol, double relax_param, double epsilon, int max_iter) {

	int sys_mat_main_diag = 4;

	if (method == JACOBI_RELAX_PARAM) {
		double* mult = mult_mat_diag_vec(sys_mat, sys_sol, n, m);

		for (int i = 0; i < n; i++) {
			sys_sol[i] = sys_sol[i] + (relax_param / sys_mat[i][sys_mat_main_diag] * (sys_right_vec[i] - mult[i]));
		}
	}
	else {

		double mult;

		for (int i = 0; i < n; i++) {
			mult = scal_mult_diag(sys_mat[i], i, sys_sol, n, m);

			sys_sol[i] = sys_sol[i] + (relax_param / sys_mat[i][sys_mat_main_diag] * (sys_right_vec[i] - mult));
		}
	}
}

double* solve_sys(METHOD method, double** sys_mat, int n, int m, double* sys_right_vec, double* sys_sol_begin_aprox, int& total_iter, double relax_param, double epsilon, int max_iter) {
	double* sys_sol = sys_sol_begin_aprox;

	total_iter = 0;
	double rel;
	do {
		iter_step(method, sys_mat, n, m, sys_right_vec, sys_sol, relax_param, epsilon, max_iter);

		rel = relative_discrepancy(sys_mat, n, m, sys_right_vec, sys_sol);

		total_iter++;
	} while (rel > epsilon and total_iter < max_iter);

	return sys_sol;
}


double* mult_mat_LDU_bandform_vec(double*** mat_elems, double* vec_elems, int n) {
	double* mult_elems = new double[n] {};

	int half_width = 1;
	int L = 0;
	int D = 1;
	int U = 2;

	for (int i = 0; i < n; i++) {
		for (int j = i - half_width; j <= i; j++) {
			if (j >= 0) {
				if (j == i)
					mult_elems[i] += mat_elems[D][0][i] * vec_elems[i];
				else {
					mult_elems[i] += mat_elems[L][i][j - (i - half_width)] * vec_elems[j];
					mult_elems[j] += mat_elems[U][i][j - (i - half_width)] * vec_elems[i];
				}
			}
		}
	}

	return mult_elems;
}

void make_LDU(double** sys_mat, int dim, int halfs_width) {
	double sum = 0;
	double mult = 0;

	for (int i = 0; i < dim; i++) {
		for (int k = 0; k < i; k++) {
			sys_mat[i][i] -= sys_mat[i][k] * sys_mat[k][k] * sys_mat[k][i];
		}

		if (i + 1 <= dim - 1) {
			for (int j = 0; j < i + 1; j++) {
				for (int k = 0; k < j; k++) {
					mult = sys_mat[i + 1][k] * sys_mat[k][k] * sys_mat[k][i + 1];
					sys_mat[j][i + 1] -= mult;
					sys_mat[i + 1][j] -= mult;
				}
				sys_mat[j][i + 1] /= sys_mat[j][j];

				sys_mat[i + 1][j] /= sys_mat[j][j];
			}
		}
		sum = 0;
		mult = 0;
	}
}

void solve_L(double** sys_mat_elems, double* sys_right_side_vec_elems, int n) {
	int half_width = 1;

	int L = 1;
	for (int i = 0; i < n; i++) {
		for (int k = 0; k < i; k++) {
			sys_right_side_vec_elems[i] -= sys_mat_elems[i][k] * sys_right_side_vec_elems[k];
		}
	}
}

void solve_D(double** sys_mat_elems, double* sys_right_side_vec_elems, int n) {
	int half_width = 1;

	int D = 0;
	for (int i = 0; i < n; i++) {
		sys_right_side_vec_elems[i] /= sys_mat_elems[i][i];
	}
}

void solve_U(double** sys_mat_elems, double* sys_right_side_vec_elems, int n) {
	int half_width = 1;

	for (int i = n - 1; i >= 0; i--) {
		for (int k = i + 1; k < n; k++) {
			sys_right_side_vec_elems[i] -= sys_mat_elems[i][k] * sys_right_side_vec_elems[k];
		}
	}
}

double* solve_sys_LDU(double** sys_mat_elems, double* sys_right_side_vec_elems, int n) {
	solve_L(sys_mat_elems, sys_right_side_vec_elems, n);
	solve_D(sys_mat_elems, sys_right_side_vec_elems, n);
	solve_U(sys_mat_elems, sys_right_side_vec_elems, n);

	return sys_right_side_vec_elems;
}

double** solve_block_relax(double**** sys_mat, int n, int block_dim, int diag_block_half_width, double** sys_right_vec, double** sys_sol_begin_aprox, double relax_param, int& total_iter, double epsilon, int max_iter) {
	double** sys_sol = sys_sol_begin_aprox;

	int sys_sol_blocks_amount = n / block_dim;

	int block_sys_mat_dim = n / block_dim;
	int half_width = 1;

	double** mult_diag_blocks = new double* [block_sys_mat_dim];

	for (int block_i = 0, diag_block = 0; block_i < block_sys_mat_dim; block_i++) {

		make_LDU(sys_mat[block_i][block_i], block_dim, diag_block_half_width);
	}

	double* mult;

	total_iter = 0;
	double rel;
	do {
		for (int block_i = 0; block_i < block_sys_mat_dim; block_i++) {

			mult = new double[block_dim] {};
			mult = vec_sum(mult, mult_mat_LDU_vec(sys_mat[block_i][block_i], sys_sol[block_i], block_dim), block_dim);

			for (int block_j = 0; block_j < block_sys_mat_dim; block_j++) {
				if (block_i != block_j)
					mult = vec_sum(mult, mult_mat_vec(sys_mat[block_i][block_j], sys_sol[block_j], block_dim), block_dim);
			}

			sys_sol[block_i] = vec_sum(solve_sys_LDU(sys_mat[block_i][block_i],
				mult_num_vec(relax_param, vec_substr(sys_right_vec[block_i],
					mult, block_dim), block_dim), block_dim), sys_sol[block_i], block_dim);

		}

		rel = block_relative_discrepancy(sys_mat, n, block_dim, sys_right_vec, sys_sol);
		/*std::cout << "iter: " << total_iter << "; rel: " << rel << std::endl;*/

		total_iter++;
	} while (rel > epsilon and total_iter <= max_iter);

	return sys_sol;
}

void cond_research_diag(METHOD diag_method) {
	std::string sys_params_file = "sys_params.txt";
	std::string sys_mat_file = "sys_mat.txt";
	std::string sys_right_vec_file = "sys_right_vec.txt";
	std::string sys_sol_begin_aprox_file = "sys_sol_begin_aprox.txt";

	double** sys_mat;
	double* sys_right_vec;
	double* sys_sol_begin_aprox;

	double epsilon;
	int n, m, max_iter;

	input_diag(sys_mat, sys_right_vec, sys_sol_begin_aprox, n, m, epsilon, max_iter, sys_params_file, sys_mat_file, sys_right_vec_file, sys_sol_begin_aprox_file);

	double step = 1e-2;
	double relax_param = 0;
	double* sys_sol = NULL;
	double* analit_sys_sol = new double[n];

	for (int i = 0; i < n; i++) {
		analit_sys_sol[i] = i + 1;
	}

	double* err;
	int total_iter;

	std::map<int, double> total_iter_relax_param;

	std::ofstream out("cond_research.txt");


	for (int i = 1; true; i++) {

		relax_param = step * i;

		sys_sol = solve_sys(diag_method, sys_mat, n, m, sys_right_vec, sys_sol_begin_aprox, total_iter, relax_param, epsilon, max_iter);

		double rel = relative_discrepancy(sys_mat, n, m, sys_right_vec, sys_sol);

		if (std::isnan(rel) or std::isinf(rel))
			break;

		err = vec_substr(analit_sys_sol, sys_sol, n);

		out << std::fixed << std::setprecision(2);

		out << relax_param << "|";

		out << std::setprecision(14) << std::scientific;

		out << sys_sol[0];

		for (int i = 1; i < n; i++) {
			out << ", " << sys_sol[i];
		}

		out << "|" << err[0];

		for (int i = 1; i < n; i++) {
			out << ", " << err[i];
		}

		out << "|" << total_iter << std::endl;

		total_iter_relax_param.insert({ total_iter,relax_param });

		delete[] sys_sol;

		sys_sol_begin_aprox = new double[n] {};
	}

	std::cout << "relax param: " << total_iter_relax_param.begin()->second << " " << total_iter_relax_param.begin()->first << std::endl;

	delete[] sys_sol;
	sys_sol_begin_aprox = new double[n] {};

	sys_sol = solve_sys(diag_method, sys_mat, n, m, sys_right_vec, sys_sol_begin_aprox, total_iter, 0.79, epsilon, max_iter);

	double rel1 = vec_norm(vec_substr(sys_sol, analit_sys_sol, n), n) / vec_norm(analit_sys_sol, n);

	double rel2 = relative_discrepancy(sys_mat, n, m, sys_right_vec, sys_sol);

	std::cout << std::setprecision(13) << std::scientific;

	std::cout << "cond num: " << rel1 / rel2;


}

void cond_research_block() {
	std::string sys_params_file = "sys_params.txt";
	std::string sys_mat_file = "sys_mat.txt";
	std::string sys_right_vec_file = "sys_right_vec.txt";
	std::string sys_sol_begin_aprox_file = "sys_sol_begin_aprox.txt";

	double**** sys_mat;
	double** sys_right_vec;
	double** sys_sol_begin_aprox;

	double epsilon;
	int n, max_iter;

	int block_dim = 4;

	input_block(sys_mat, n, sys_right_vec, sys_sol_begin_aprox, epsilon, max_iter, block_dim, sys_params_file, sys_mat_file, sys_right_vec_file, sys_sol_begin_aprox_file);

	int block_mat_dim = n / block_dim;

	double step = 1e-2;
	double relax_param = 0;
	double** sys_sol;
	double** analit_sys_sol = new double* [block_mat_dim];

	int num = 1;
	for (int block_i = 0; block_i < block_mat_dim; block_i++) {
		analit_sys_sol[block_i] = new double[block_dim];

		for (int i = 0; i < block_dim; i++, num++) {
			analit_sys_sol[block_i][i] = num;
		}
	}

	double* err = new double[n];
	int total_iter;

	std::map<int, double> total_iter_relax_param;

	std::ofstream out("cond_research.txt");


	for (int i = 1; true; i++) {

		relax_param = step * i;

		sys_sol = solve_block_relax(sys_mat, n, block_dim, 0, sys_right_vec, sys_sol_begin_aprox, relax_param, total_iter, epsilon, max_iter);

		double rel = block_relative_discrepancy(sys_mat, n, block_dim, sys_right_vec, sys_sol);

		if (std::isnan(rel) or std::isinf(rel))
			break;

		int err_i = 0;
		for (int block_i = 0; block_i < block_mat_dim; block_i++) {
			for (int i = 0; i < block_dim; i++, err_i++) {
				err[err_i] = analit_sys_sol[block_i][i] - sys_sol[block_i][i];
			}
		}

		out << std::fixed << std::setprecision(2);

		out << relax_param << "|";

		out << std::setprecision(14) << std::scientific;

		for (int block_i = 0; block_i < block_mat_dim; block_i++) {
			for (int i = 0; i < block_dim; i++) {
				out << " " << sys_sol[block_i][i];
			}
		}

		out << "|" << err[0];

		for (int i = 1; i < n; i++) {
			out << " " << err[i];
		}

		out << "|" << total_iter << std::endl;

		total_iter_relax_param.insert({ total_iter,relax_param });

		delete[] sys_sol;

		input_block(sys_mat, n, sys_right_vec, sys_sol_begin_aprox, epsilon, max_iter, block_dim, sys_params_file, sys_mat_file, sys_right_vec_file, sys_sol_begin_aprox_file);
	}

	std::cout << "relax param: " << total_iter_relax_param.begin()->second << " " << total_iter_relax_param.begin()->first << std::endl;

	delete[] sys_sol;

	input_block(sys_mat, n, sys_right_vec, sys_sol_begin_aprox, epsilon, max_iter, block_dim, sys_params_file, sys_mat_file, sys_right_vec_file, sys_sol_begin_aprox_file);

	sys_sol = solve_block_relax(sys_mat, n, block_dim, 0, sys_right_vec, sys_sol_begin_aprox, total_iter_relax_param.begin()->second, total_iter, epsilon, max_iter);

	double* substr;
	double sum1 = 0;
	double sum2 = 0;
	for (int block_i = 0; block_i < block_mat_dim; block_i++) {
		substr = vec_substr(sys_sol[block_i], analit_sys_sol[block_i], block_dim);
		sum1 += scal_mult(substr, substr, block_dim);
		sum2 += scal_mult(analit_sys_sol[block_i], analit_sys_sol[block_i], block_dim);
	}

	double rel1 = sqrt(sum1) / sqrt(sum2);

	double rel2 = block_relative_discrepancy(sys_mat, n, block_dim, sys_right_vec, sys_sol);

	std::cout << std::setprecision(13) << std::scientific;

	std::cout << "cond num: " << rel1 / rel2;

}

int main() {
	std::string sys_params_file = "sys_params.txt";
	std::string sys_mat_file = "sys_mat.txt";
	std::string sys_right_vec_file = "sys_right_vec.txt";
	std::string sys_sol_begin_aprox_file = "sys_sol_begin_aprox.txt";

	double**** sys_mat;
	double** sys_right_vec;
	double** sys_sol_begin_aprox;

	double epsilon;
	int n, max_iter;

	int total_iter;

	int block_dim = 3;

	input_block(sys_mat, n, sys_right_vec, sys_sol_begin_aprox, epsilon, max_iter, block_dim, sys_params_file, sys_mat_file, sys_right_vec_file, sys_sol_begin_aprox_file);

	double** sys_sol = solve_block_relax(sys_mat, n, block_dim, 0, sys_right_vec, sys_sol_begin_aprox, 1.5, total_iter, epsilon, max_iter);

	output_block_vec(sys_sol, n, block_dim, "block_res.txt");
}


