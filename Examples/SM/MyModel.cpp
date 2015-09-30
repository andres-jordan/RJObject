#include "MyModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include <cmath>

using namespace std;
using namespace Eigen;
using namespace DNest3;

MyModel::MyModel()
:objects(3, 10, false, MyDistribution(-6., 1., 1E-3, 1E3, -6, 6))
,mu(Data::get_instance().get_t().size())
,C(Data::get_instance().get_t().size(),
       Data::get_instance().get_t().size())
{

}

void MyModel::fromPrior()
{
	objects.fromPrior();
	objects.consolidate_diff();
	sigma = exp(log(1E-3) + log(1E6)*randomU());
	calculate_mu();
	calculate_C();
}

void MyModel::calculate_C()
{
        // Get the data
        const vector<double>& t = Data::get_instance().get_t();
        const vector<double>& sig = Data::get_instance().get_sig();

        // Update or from scratch?
        bool update = (objects.get_added().size() < objects.get_components().size()) &&
                        (staleness <= 10);


 		// Get the components
		const vector< vector<double> >& components = (update)?(objects.get_added()):
					(objects.get_components());


		double A, s, v;
        for(size_t i=0; i<Data::get_instance().get_t().size(); i++)
        {
                for(size_t j=i; j<Data::get_instance().get_t().size(); j++)
                {
                		C(i,j) = 0;
                		for(size_t m=1; m<components.size(); m++){
                			A = exp(components[m][0]);
                			v = exp(components[m][1]);
                            s = exp(components[m][2]);
                			C(i,j) += (A*cos(2*M_PI*s*(t[i]-t[j])) * exp(-2*M_PI*M_PI*(t[i]-t[j])*(t[i]-t[j])*v));
                		}           
                        if(i==j)
                                C(i, j) += sig[i]*sig[i] + sigma*sigma;
                        else
                                C(j, i) = C(i, j);
                }
        }
}


void MyModel::calculate_mu()
{
	// Get the times from the data
	const vector<double>& t = Data::get_instance().get_t();

	// Update or from scratch?
	bool update = false;//(objects.get_added().size() < objects.get_components().size());

	// Zero the signal
	if(!update)
		mu.assign(mu.size(), 0.);

	// at the moment we do *not* want anything at mu
	//mu[i] += A*sin(2.*M_PI*t[i]/T + phi);
	for(size_t i=0; i<t.size(); i++){
		mu[i] += 0;
	}

	
}

double MyModel::perturb()
{
	double logH = 0.;

	if(randomU() <= 0.75)
	{
		logH += objects.perturb();
		objects.consolidate_diff();
		calculate_mu();
	}
	else
	{
		sigma = log(sigma);
		sigma += log(1E6)*randh();
		sigma = mod(sigma - log(1E-5), log(1E6)) + log(1E-5);
		sigma = exp(sigma);
	}
	calculate_C();

	return logH;
}

double MyModel::logLikelihood() const
{
        // Get the data
        const vector<double>& y = Data::get_instance().get_y();

        VectorXd residual(y.size());
        for(size_t i=0; i<y.size(); i++)
                residual(i) = y[i] - mu[i];

        Eigen::LLT<Eigen::MatrixXd> cholesky = C.llt();
        MatrixXd L = cholesky.matrixL();

        double logDeterminant = 0.;
        for(size_t i=0; i<y.size(); i++)
                logDeterminant += 2.*log(L(i,i));

        // C^-1*(y-mu)
        VectorXd solution = cholesky.solve(residual);

        // y*solution
        double exponent = 0.;
        for(size_t i=0; i<y.size(); i++)
                exponent += residual(i)*solution(i);

        double logL = -0.5*y.size()*log(2*M_PI)
                                        - 0.5*logDeterminant - 0.5*exponent;

        if(isnan(logL) || isinf(logL))
                logL = -1E300;

        return logL;
}


void MyModel::print(std::ostream& out) const
{
	for(size_t i=0; i<mu.size(); i++)
		out<<mu[i]<<' ';
	out<<sigma<<' ';
	objects.print(out); out<<' ';
}

string MyModel::description() const
{
	return string("objects, sigma");
}

