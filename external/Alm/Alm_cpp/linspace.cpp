#include <Eigen/Dense>
#include <iostream>
#include <vector>

using Eigen::VectorXd;


// Taken from https://stackoverflow.com/questions/27028226/python-linspace-in-c
std::vector<double> linspace_vec(const long double start, const long double end, const long num)
{
  std::vector<double> linspaced;

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}
// -----

// My linspace function
VectorXd linspace(const long double start_in, const long double end_in, const long num_in)
{
	if (num_in == 0) {
		std::cout << " num_in in linspace is 0. Cannot create a linspace vector. Returning -1." << std::endl;
		VectorXd linspaced(1);
		linspaced[0]=-1;
		return linspaced;
	}
	VectorXd linspaced(num_in);

	const long double delta = (end_in - start_in) / (num_in - 1);
	for(long i=0 ; i< num_in ; i++){
		linspaced[i]=start_in + delta*i;
	}

	/*std::cout << "start_in =" << start_in << std::endl;
	std::cout << "end_in =" << end_in << std::endl;
	std::cout << "num_in =" << num_in << std::endl;
	std::cout << "delta =" << delta << std::endl;
	std::cout << "linspaced.size() =" << linspaced.size() << std::endl;
	
	for(long i=0; i<linspaced.size(); i++)
	{
		std::cout << linspaced[i] << std::endl;
	}
	std::cout << "linspace needs a thorough test" << std::endl;
	std::cout << "exiting now"<< std::endl;
	exit(EXIT_SUCCESS);
	*/
	return linspaced;
}
