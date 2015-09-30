#include "MyDistribution.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace DNest3;

MyDistribution::MyDistribution(double a_min, double a_max,
					double mu_min, double mu_max, double v_min, double v_max)
:a_min(a_min)
,a_max(a_max)
,mu_min(mu_min)
,mu_max(mu_max)
,v_min(v_min)
,v_max(v_max)
{

}

void MyDistribution::fromPrior()
{
	mu = exp(log(mu_min) + log(mu_max/mu_min)*randomU());
}

double MyDistribution::perturb_parameters()
{
	double logH = 0.;

	mu = log(mu);
	mu += log(mu_max/mu_min)*pow(10., 1.5 - 6.*randomU())*randn();
	mu = mod(mu - log(mu_min), log(mu_max/mu_min)) + log(mu_min);
	mu = exp(mu);

	return logH;
}

// vec[0] = "position" (log-period) --> amplitude
// vec[1] = amplitude --> inverse variance 
// vec[2] = phase     --> "position" (log frequency)

double MyDistribution::log_pdf(const std::vector<double>& vec) const
{
	if(vec[0] < a_min || vec[0] > a_max || vec[1] < 0. ||
			vec[2] < v_min || vec[2] > v_max)
		return -1E300;

	return -log(mu) - vec[1]/mu;
}

void MyDistribution::from_uniform(std::vector<double>& vec) const
{
	vec[0] = a_min + (a_max - a_min)*vec[0];
	vec[1] = -mu*log(1. - vec[1]);
	vec[2] = v_min + (v_max - v_min)*vec[2];
}

void MyDistribution::to_uniform(std::vector<double>& vec) const
{
	vec[0] = (vec[0] - a_min)/(a_max - a_min);
	vec[1] = 1. - exp(-vec[1]/mu);
	vec[2] = (vec[2] - v_min)/(v_max - v_min);
}

void MyDistribution::print(std::ostream& out) const
{
	out<<mu<<' ';
}

