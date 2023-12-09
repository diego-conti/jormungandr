#include "labeled_tree.h"
#include "weightmatrix.h"
#include "diagramprocessor.h"
#include "nicediagramsinpartition.h"
#include "negativeindices.h"
#include "injectivemaps.h"
#include "niceliegroupwithmetric.h"
#include "polynomialsolver.h"
#include "restrictionsandoptions.h"
#include "generalizednilsoliton.h"
#include "nondiagonalmatrix.h"
#include "nondiagonalgeneralizednilsolitons.h"
#include "splitter.h"
#include "partitioned.h"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

struct surjective_type_t {} surjective_type;

string diagram_name(const LabeledTree& diagram) {
	stringstream s;
	s<<horizontal(lower_central_series(diagram),"")<<":"<<diagram.number();
	return s.str();
}

void study_diagram(LabeledTree diagram, Restrictions restrictions, OutputOptions options, surjective_type_t) {		
	NondiagonalGeneralizedNilsolitons n{diagram,restrictions,options};	
	auto section=n.section_of_generalized_nilsolitons_surjective_type();
	auto partitioned = PartitionedGeneralizedNilsolitonSurjectiveType{section};	
	partitioned.print_nilsolitons(cout,diagram_name(diagram));
	partitioned.print_einstein_solvmanifolds(*options.metrics_stream);
}

void study_diagram(LabeledTree diagram, Restrictions restrictions, OutputOptions options) {		
	NondiagonalGeneralizedNilsolitons N{diagram,restrictions,options};
	cout<<"diagram: "<<diagram.as_string()<<endl;
	for (auto n: N.generalized_nilsolitons()) 
		n.print(cout);	
}

unique_ptr<LabeledTree> read_diagram(const string& s) {
	stringstream str{s};
	return LabeledTree::from_stream(str);
}

void print_help(string program_name, const po::options_description& desc) {
	cout<<"Jǫrmungandr"<<endl;
	cout<<"A program to costruct Einstein solvmanifolds that are not of Iwasawa type"<<endl;
	cout << "Usage: "<< program_name<< " [options]\n";
	cout << desc;
}

struct PartitionRange {
	vector<int> first;
	vector<int> last;
	list<vector<int>> partitions_in_range(int dimension) const {
		auto all_partitions=partitions(dimension);
		auto i=all_partitions.begin();
		if (!first.empty()) 
			while (*i!=first) ++i;
		auto j=i;
		if (last.empty()) j=all_partitions.end();
		else
			while ((*j++)!=last) ++j;
		return {i,j};
	}
};


class Runner {
	const po::variables_map& command_line_variables;
	Restrictions restrictions;
	OutputOptions options;
	string metrics_filename;
	PartitionRange partition_range;
	void print_metrics_filename() const {
		if (!metrics_filename.empty()) cout<<"metrics written to "<<metrics_filename<<endl;	
	}
	template<typename... Tags>
	int enumerate_nice_diagrams_in_partition(vector<int> partition, Tags... tags) {
		for (auto diagram: nice_diagrams_in_partition(partition)) {
			diagram.invert_nodes();		
			study_diagram(diagram,restrictions, options, tags...);
		}
		return 0;	
	}
	template<typename... Tags>
	void enumerate_nice_diagrams(int dimension, Tags... tags) {
		for (auto& partition: partition_range.partitions_in_range(dimension)) {
			enumerate_nice_diagrams_in_partition(partition,tags...);
		}
	}
	int enumerate_nice_diagrams_and_print(int dimension, surjective_type_t) {
		cout<<"\\begin{array}{cccc}"<<endl<<"\\Delta&\\lie g & D & S\\\\"<<endl;
		enumerate_nice_diagrams(dimension,surjective_type);
		cout<<"\\end{array}"<<endl;
		print_metrics_filename();
		return 0;
	}
	int enumerate_nice_diagrams_and_print(int dimension) {
		enumerate_nice_diagrams(dimension);		
		return 0;
	}
	int parse_and_study_diagram(const string& d) {
		try {
			auto diagram=read_diagram(d);
			if (diagram) {
				if (diagram->weight_basis({}).properties().dimension_cokernel_M_Delta()==0) {
					study_diagram(*diagram,restrictions,options,surjective_type);
					print_metrics_filename();
				}
				else
					 study_diagram(*diagram,restrictions,options);			
				return 0;
			}
			else {
				cerr<<"cannot parse diagram"<<endl;
				return 1;		
			}
		}

		catch (std::exception& e) {
			cerr<<e.what()<<endl;			
			return 1;	
		}				
	}
	bool verbos=false;
public:
	Runner(const po::variables_map& command_line_variables) : command_line_variables{command_line_variables} {
		if (command_line_variables.count("free-parameters")) restrictions=restrictions.free_parameters(command_line_variables["free-parameters"].as<string>());	
		if (command_line_variables.count("offdiagonal-elements")) restrictions=restrictions.bound_on_A(command_line_variables["offdiagonal-elements"].as<string>());	
		if (command_line_variables.count("verbose")) options.verbose=true;
		if (command_line_variables.count("output-metrics")) {
			if (!command_line_variables.count("surjective"))
				cerr<<"WARNING: --output-metrics only effective if --surjective is indicated"<<endl;
			else {
				metrics_filename=command_line_variables["output-metrics"].as<string>();
				options.metrics_stream=make_shared<ofstream>(metrics_filename);			
			}
		}
	}
	int run() {
		bool surjective=command_line_variables.count("surjective");
		if (command_line_variables.count("from-partition"))
			partition_range.first=command_line_variables["from-partition"].as<vector<int>>();
		if (command_line_variables.count("to-partition"))
			partition_range.last=command_line_variables["to-partition"].as<vector<int>>();

		if (command_line_variables.count("dimension"))
			return surjective?
				 enumerate_nice_diagrams_and_print(command_line_variables["dimension"].as<int>(),surjective_type)
				:			
				 enumerate_nice_diagrams_and_print(command_line_variables["dimension"].as<int>());
		if (command_line_variables.count("partition"))
			return surjective?
				 enumerate_nice_diagrams_in_partition(command_line_variables["partition"].as<vector<int>>(),surjective_type)
				:			
				 enumerate_nice_diagrams_in_partition(command_line_variables["partition"].as<vector<int>>());				 
		else if (command_line_variables.count("digraph")) 
			return parse_and_study_diagram(command_line_variables["digraph"].as<string>());
		else return 1;		
	}
};


int main(int argc, char** argv) {	
	cout<<latex;	
 	po::options_description desc("Allowed options");
    desc.add_options()
		("help", "produce help message")
        ("dimension", po::value<int>(),"process all nice diagrams of dimension <arg>")
        ("partition", po::value<vector<int>>()->multitoken(),"process all nice diagrams associated to the partition <arg>")
        ("from-partition", po::value<vector<int>>()->multitoken(),"process partitions of a fixed dimension starting from <arg>")
        ("to-partition", po::value<vector<int>>()->multitoken(),"process partitions of a fixed dimensions until <arg>")
        ("digraph", po::value<string>(),
              "only process diagram indicated by <arg>")
		("verbose", "print extra output")
        ("surjective", "only study situations in which the root matrix is surjective")
		("output-metrics", po::value<string>(), "output the list of Einstein metrics obtained to <arg>")
		("free-parameters", po::value<string>(), "only consider cases in which the remaining number of free parameters A_ij satisfies <arg>, where <arg> takes the form =n, <n or >n")
		("offdiagonal-elements", po::value<string>(), "only consider sets of indices A of size satisfying <arg>, where <arg> takes the form =n, <n or >n");	
	
 	po::variables_map vm;
	try {
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
	}
	catch (const boost::program_options::unknown_option& e) {
		print_help(argv[0],desc);
		return 1;
	}
	if (vm.count("help") || Runner{vm}.run()) print_help(argv[0],desc);	
}
