#include "ExampleFunction.h"
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#define weight_eta 200

ExampleFunction::ExampleFunction(Placement &placement) : _placement(placement)
{
    NumModules = _placement.numModules();

    double width = _placement.boundryRight() - _placement.boundryLeft();
    double height = _placement.boundryTop() - _placement.boundryBottom();
	double planArea = width * height;

	double totalArea = 0;
	for (int id = 0; id < NumModules; id++)
		totalArea += _placement.module(id).area();

	Tb = totalArea / planArea;

	double boundryLeft = _placement.boundryLeft();
	double boundryBottom = _placement.boundryBottom();

	wb = width / NumBinPerRow;
	hb = height / NumBinPerRow;

	for (int i = 0; i < NumBinPerRow; i++) {
		double yi = boundryBottom + (0.5 + i) * hb;
		for (int j = 0; j < NumBinPerRow; j++) {
			double xi = boundryLeft + (0.5 + j) * wb;
			bin_matrix[i][j].xi = xi;
			bin_matrix[i][j].yi = yi;
		}
	}

    eta = (width + height) / weight_eta;
	//cout << "eta:" << eta << endl;
    Ci = 0;
}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{
    g = vector<double>(g.size(), 0); //initialize g


	double LSE = 0;

	vector<double> _exp(NumModules * 4);
	for (int i = 0; i < NumModules; i++) {
		_exp[i * 4] = exp(x[i * 2] / eta);
		_exp[i * 4 + 1] = exp(-x[i * 2] / eta);
		_exp[i * 4 + 2] = exp(x[i * 2 + 1] / eta);
		_exp[i * 4 + 3] = exp(-x[i * 2 + 1] / eta);
	}

	int NumNets = _placement.numNets();
	for (int i_net = 0; i_net < NumNets; i_net++) {
		Net &net = _placement.net(i_net);

		double LSE_x = 0;
		double LSE_mx = 0;
		double LSE_y = 0;
		double LSE_my = 0;
		int NumPins = net.numPins();
		for (int i_pin = 0; i_pin < NumPins; i_pin++) {
			int id = net.pin(i_pin).moduleId();
			LSE_x += _exp[id * 4];
			LSE_mx += _exp[id * 4 + 1];
			LSE_y += _exp[id * 4 + 2];
			LSE_my += _exp[id * 4 + 3];
		}
		//log sum
		LSE += log(LSE_x) + log(LSE_mx) + log(LSE_y) + log(LSE_my);
		//cout << "LSE*eta " << eta * LSE;

		//grad update
		for (int i_pin = 0; i_pin < NumPins; i_pin++) {
			int id = net.pin(i_pin).moduleId();
			g[id * 2] += _exp[id * 4] / LSE_x;
			g[id * 2] -= _exp[id * 4 + 1] / LSE_mx;
			g[id * 2 + 1] += _exp[id * 4 + 2] / LSE_y;
			g[id * 2 + 1] -= _exp[id * 4 + 3] / LSE_my;
		}
	}

	double Db = 0;

	for (int i = 0; i < NumBinPerRow; i++) {
		for (int j = 0; j < NumBinPerRow; j++) {
			vector<double> _g(g.size(), 0);

			double bin_density = 0;
			for (int id = 0; id < NumModules; id++) {
				double wi = _placement.module(id).width();
				double hi = _placement.module(id).height();
				//Db on x
				double Dbx = 0;
				double Dbx_grad = 0;
				double dx = x[id * 2] - bin_matrix[i][j].xi;
				double dx_abs = abs(dx);
				if (dx_abs <= wb / 2 + wi / 2) {
					double a = 4 / ((wb + wi) * (2 * wb + wi));
					Dbx = 1 - a * dx_abs * dx_abs;
					Dbx_grad = -2 * a * dx;
				}
				else if (dx_abs <= wb + wi / 2) {
					double b = 4 / (wb * (2 * wb + wi));
					Dbx = b * (dx_abs - wb - wi / 2) * (dx_abs - wb - wi / 2);
					if (dx > 0)
						Dbx_grad = 2 * b * (dx - wb - wi / 2) * 1;
					else
						Dbx_grad = 2 * b * (dx - wb - wi / 2) * -1;
				}
				else {
					Dbx = 0;
					Dbx_grad = 0;
				}
				//Db on y
				double Dby = 0;
				double Dby_grad = 0;
				double dy = x[id * 2 + 1] - bin_matrix[i][j].yi;
				double dy_abs = abs(dy);
				if (dy_abs <= hb / 2 + hi / 2) {
					double a = 4 / ((hb + hi) * (2 * hb + hi));
					Dby = 1 - a * dy_abs * dy_abs;
					Dby_grad = -2 * a * dy;
				}
				else if (dy_abs <= hb + hi / 2) {
					double b = 4 / (hb * (2 * hb + hi));
					Dby = b * (dy_abs - hb - hi / 2) * (dy_abs - hb - hi / 2);
					if (dy > 0)
						Dby_grad = 2 * b * (dy - hb - hi / 2) * 1;
					else
						Dby_grad = 2 * b * (dy - hb - hi / 2) * -1;
				}
				else {
					Dby = 0;
					Dby_grad = 0;
				}

				bin_density += Dbx * Dby;

				_g[id * 2] = Dbx_grad * Dby;
				_g[id * 2 + 1] = Dbx * Dby_grad;
			}

			double delta_Db = max(bin_density - Tb, 0.0);
			Db += delta_Db * delta_Db;

			//grad update
			for (int id = 0; id < NumModules; id++) {
				g[id * 2] += Ci * 2 * delta_Db * _g[id * 2];
				g[id * 2 + 1] += Ci * 2 * delta_Db * _g[id * 2 + 1];
			}
		}
	}


	f = eta * LSE + Ci * Db;
}

void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{	//cal LSE
	double LSE = 0;

	vector<double> _exp(NumModules * 4);
	for (int i = 0; i < NumModules; i++) {
		_exp[i * 4] = exp(x[i * 2] / eta);
		_exp[i * 4 + 1] = exp(-x[i * 2] / eta);
		_exp[i * 4 + 2] = exp(x[i * 2 + 1] / eta);
		_exp[i * 4 + 3] = exp(-x[i * 2 + 1] / eta);
	}

	int NumNets = _placement.numNets();
	for (int i_net = 0; i_net < NumNets; i_net++) {
		Net &net = _placement.net(i_net);

		double LSE_x = 0;
		double LSE_mx = 0;
		double LSE_y = 0;
		double LSE_my = 0;
		int NumPins = net.numPins();
		for (int i_pin = 0; i_pin < NumPins; i_pin++) {
			int id = net.pin(i_pin).moduleId();

			LSE_x += _exp[id * 4];
			LSE_mx += _exp[id * 4 + 1];
			LSE_y += _exp[id * 4 + 2];
			LSE_my += _exp[id * 4 + 3];
		}

		LSE += log(LSE_x) + log(LSE_mx) + log(LSE_y) + log(LSE_my);
	}

	//return eta * LSE;

    //double lse = eta * LSE;

	//cal Db
	double Db = 0;

	for (int i = 0; i < NumBinPerRow; i++) {
		for (int j = 0; j < NumBinPerRow; j++) {
			double bin_density = 0;
			for (int id = 0; id < NumModules; id++) {
				double wi = _placement.module(id).width();
				double hi = _placement.module(id).height();
				//Db on x
				double Dbx = 0;
				double dx = abs(x[id * 2] - bin_matrix[i][j].xi);
				if (dx <= wb / 2 + wi / 2) {
					double a = 4 / ((wb + wi) * (2 * wb + wi));
					Dbx = 1 - a * dx * dx;
				}
				else if (dx <= wb + wi / 2) {
					double b = 4 / (wb * (2 * wb + wi));
					Dbx = b * (dx - wb - wi / 2) * (dx - wb - wi / 2);
				}
				else {
					Dbx = 0;
				}


				//Db on y
				double Dby = 0;
				double dy = abs(x[id * 2 + 1] - bin_matrix[i][j].yi);
				if (dy <= hb / 2 + hi / 2) {
					double a = 4 / ((hb + hi) * (2 * hb + hi));
					Dby = 1 - a * dy * dy;
				}
				else if (dy <= hb + hi / 2) {
					double b = 4 / (hb * (2 * hb + hi));
					Dby = b * (dy - hb - hi / 2) * (dy - hb - hi / 2);
				}
				else {
					Dby = 0;
				}

				bin_density += Dbx * Dby;
			}

			double delta_Db = max(bin_density - Tb, 0.0);
			Db += delta_Db * delta_Db;
		}
	}

    f = eta * LSE + Ci * Db;
}

unsigned ExampleFunction::dimension()
{
    return 2*NumModules; 
    // each two dimension represent the X and Y dimensions of each block
}
