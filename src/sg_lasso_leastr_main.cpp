#include <argparse.hpp>
#include <armadillo>
#include "sg_lasso_leastr.hpp"

using namespace std;
using namespace arma;


map<string, string> processSlepOpts(string filename);
std::vector<std::tuple<double, double>> readLambdaList(string filename);
std::string lambdaLabel(const double arr[], int idx);
string writeModelToXMLStream(string);

int main(int argc, char *argv[]) {
  argparse::ArgumentParser program("sg_lasso_leastr");

  program.add_argument("-f", "--features")
    .required()
    .help("Specify the input features file.");

  program.add_argument("-n", "--groups")
    .required()
    .help("Specify the group indices file.");

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

  program.add_argument("-c", "--gene_count_threshold")
    .default_value(0)
    .help("Specify gene selection cutoff for grid-based model building.")
    .scan<'i', int>();

  try {
    program.parse_args(argc, argv);
  }
  catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    std::exit(1);
  }

  //auto input = program.get<int>("square");
  //std::cout << (input * input) << std::endl;

  double lambda[2] = {program.get<double>("lambda1"), program.get<double>("lambda2")};

  mat features;
  mat opts_ind;
  rowvec responses;

  SGLassoLeastR* sgl;

  int auto_cancel_threshold = program.get<int>("gene_count_threshold");

  //features.load(csv_name(program.get<std::string>("features"),csv_opts::semicolon));
  features.load(csv_name(program.get<std::string>("features"),csv_opts::trans));

  //responses.load(csv_name(program.get<std::string>("response"),csv_opts::semicolon));
  responses.load(csv_name(program.get<std::string>("response"),csv_opts::trans));





  if (responses.n_cols != features.n_cols)
  {
    //Log::Fatal << "The responses must have the same number of columns as the feature set." << endl;
    throw std::invalid_argument("\nThe responses must have the same number of columns as the feature set.\n");
  }

  opts_ind.load(csv_name(program.get<std::string>("groups"),csv_opts::trans));

  if (program.get<std::string>("lambda_list") != "-") {
    std::vector<std::tuple<double, double>> lambda_list = readLambdaList(program.get<std::string>("lambda_list"));
    double max_glambda2 = 1.0;
    double min_glambda2 = 1.0;
    for (const auto& item : lambda_list) {
	  double glambda[2] = {std::get<0>(item), std::get<1>(item)};
	  if (min_glambda2 > glambda[1]) {
        min_glambda2 = glambda[1];
      }
	  std::cout << glambda[0] << " - " << glambda[1] << std::endl;
	  if (glambda[1]>max_glambda2) {
        std::cout << "Skipping this lambda value due to gene count thresholding..." << std::endl;
        continue;
      }
	  sgl = new SGLassoLeastR(features, responses, opts_ind, glambda, processSlepOpts(program.get<std::string>("slep")));

      //TODO: make out filename reflect lambda pair
	  ofstream fileStream(program.get<std::string>("output") + lambdaLabel(glambda,0) + lambdaLabel(glambda,1) + ".xml");
	  if (fileStream.is_open())
	  {
        sgl->writeModelToXMLStream(fileStream);
        fileStream.close();
      } else {
        std::cout << "Could not open output file for writing." << std::endl;
      }
      std::cout << "Non-zero gene count: " << sgl->NonZeroGeneCount() << std::endl;
      if (sgl->NonZeroGeneCount()<=auto_cancel_threshold){
        if (glambda[1] == min_glambda2) {
          std::cout << "Skipping all further lambda pairs due to gene count thresholding..." << std::endl;
          break;
        }
        max_glambda2 = glambda[1];
      }
    }

    return 0;

  }



  sgl = new SGLassoLeastR(features, responses, opts_ind, lambda, processSlepOpts(program.get<std::string>("slep")));

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

