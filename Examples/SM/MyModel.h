#ifndef _MyModel_
#define _MyModel_

#include "Model.h"
#include <vector>
#include <RJObject.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include "MyDistribution.h"

class MyModel:public DNest3::Model
{
	private:
		RJObject<MyDistribution> objects;
		double sigma; // Noise standard deviation

		// The signal
		std::vector<long double> mu;
		void calculate_mu();

        // The covariance matrix for the data
        Eigen::MatrixXd C;
        void calculate_C();

        unsigned int staleness;

	public:
		MyModel();

		// Generate the point from the prior
		void fromPrior();

		// Metropolis-Hastings proposals
		double perturb();

		// Likelihood function
		double logLikelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;
};

#endif

