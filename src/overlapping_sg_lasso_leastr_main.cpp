#include <argparse.hpp>
#include <armadillo>
#include "overlapping_sg_lasso_leastr.hpp"

using namespace std;
using namespace arma;


map<string, string> processSlepOpts(string filename);
std::vector<std::tuple<double, double>> readLambdaList(string filename);
std::string lambdaLabel(const double arr[], int idx);
string writeModelToXMLStream(string);

int main(int argc, char *argv[]) {
  argparse::ArgumentParser program("sg_lasso");

  program.add_argument("-f", "--features")
    .required()
    .help("Specify the input features file.");

  program.add_argument("-n", "--groups")
    .required()
    .help("Specify the group indices file.");

  program.add_argument("-g", "--field")
    .required()
    .help("Feature indices file for overlapping groups.");

  program.add_argument("-r", "--response")
    .required()
    .help("Specify the response file.");

  program.add_argument("-w", "--output")
    .required()
    .help("specify the output file.");

  program.add_argument("-s", "--slep")
    .default_value(std::string("-"))
    .help("Specify a file of key/value SLEP options.");

  program.add_argument("-l", "--lambda_list")
    .default_value(std::string("-"))
    .help("Specify a file of lambda value pairs.");

  program.add_argument("-z", "--lambda1")
    .default_value(0.1)
    .help("Specify individual feature sparsity.")
    .scan<'g', double>();

  program.add_argument("-y", "--lambda2")
    .default_value(0.1)
    .help("Specify group feature sparsity.")
    .scan<'g', double>();

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  double lambda[2] = {program.get<double>("lambda1"), program.get<double>("lambda2")};

  mat features;
  mat opts_ind;
  rowvec field;
  rowvec responses;

  OLSGLassoLeastR* sgl;

  features.load(csv_name(program.get<std::string>("features"),csv_opts::trans));
  responses.load(csv_name(program.get<std::string>("response"),csv_opts::trans));


  if (responses.n_cols != features.n_cols)
  {
    //Log::Fatal << "The responses must have the same number of columns as the feature set." << endl;
    throw std::invalid_argument("\nThe responses must have the same number of columns as the feature set.\n");
  }

  opts_ind.load(csv_name(program.get<std::string>("groups"),csv_opts::trans));

  //field.load(csv_name(program.get<std::string>("field"),csv_opts::semicolon + csv_opts::trans));
  field.load(csv_name(program.get<std::string>("field")));

  if (field(field.n_cols - 1) == 0){field.resize(field.n_elem - 1);}

  if (program.get<std::string>("lambda_list") != "-") {
    std::vector<std::tuple<double, double>> lambda_list = readLambdaList(program.get<std::string>("lambda_list"));

    for (const auto& item : lambda_list) {
	  double glambda[2] = {std::get<0>(item), std::get<1>(item)};
	  std::cout << glambda[0] << " - " << glambda[1] << std::endl;
	  sgl = new OLSGLassoLeastR(features, responses, opts_ind, field, glambda, processSlepOpts(program.get<std::string>("slep")));

      //TODO: make out filename reflect lambda pair
	  ofstream fileStream(program.get<std::string>("output") + lambdaLabel(glambda,0) + lambdaLabel(glambda,1) + ".xml");
	  if (fileStream.is_open())
	  {
        sgl->writeModelToXMLStream(fileStream);
        fileStream.close();
      } else {
        std::cout << "Could not open output file for writing." << std::endl;
      }
    }

    return 0;

  }



  sgl = new OLSGLassoLeastR(features, responses, opts_ind, field, lambda, processSlepOpts(program.get<std::string>("slep")));

  //std::cout << sgl->writeModelToXMLStream();
  ofstream fileStream(program.get<std::string>("output") + ".xml");
  if (fileStream.is_open())
  {
    //fileStream << sgl->writeModelToXMLStream();
    sgl->writeModelToXMLStream(fileStream);
    fileStream.close();
  } else {
    std::cout << "Could not open output file for writing." << std::endl;
  }

  return 0;

}


map<string, string> processSlepOpts(string filename)
{
  map<string, string> slep_opts;
  std::cout << "Processing SLEP options file: " << filename << "..." << std::endl;

  string line;
  ifstream optsFile (filename);

  if (optsFile.is_open())
  {

    int splitpos;
    string opt_key;
    while (getline(optsFile, line))
    {
      splitpos = line.find("\t");
      if (splitpos != std::string::npos)
      {
        opt_key = line.substr(0, line.find("\t"));
        slep_opts[opt_key] = line.substr(line.find("\t")+1, std::string::npos);
      }
    }
    optsFile.close();
  }

  return slep_opts;
}


std::vector<std::tuple<double, double>> readLambdaList(string filename)
{
  std::vector<std::tuple<double, double>> data;
  std::string line;
  std::ifstream file(filename);

  if (!file) {
    std::cerr << "Unable to open the file." << std::endl;
    return data;
  }



  while (std::getline(file, line)) {
    std::istringstream iss(line);
	float val1, val2;

	if (iss >> val1 && iss >> val2) {
	  data.push_back(std::make_tuple(val1, val2));
    }
  }

  return data;
}

std::string lambdaLabel(const double arr[], int idx) {


    // Convert the double to a string
    std::string doubleStr = std::to_string(arr[idx]);

    // Find "0." at the beginning of the string and erase it if present
    if(doubleStr.substr(0, 2) == "0.") {
        doubleStr.erase(0, 2);
    }

    size_t lastNonZero = doubleStr.rfind('0', doubleStr.size()-1);
    if (lastNonZero == std::string::npos) {  // no zeroes found
        return "_" + doubleStr;
    }

    size_t pos = lastNonZero;
    while (pos > 0 && doubleStr[pos-1] == '0') {
        --pos;
    }
    return "_" + doubleStr.substr(0, pos);
}