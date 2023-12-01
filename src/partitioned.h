
struct HashGeneralizedNilsolitonSurjectiveTypeNoSignature {
	matrix D;
	map<tuple<int,int,int>,ex> c_ijk;	//the structure constants, i.e. [e_i,e_j]=c_ij e_k, where i<j 		
	bool operator<(const HashGeneralizedNilsolitonSurjectiveTypeNoSignature& other) const {
		return ex_is_less()(D,other.D) || (D==other.D && c_ijk<other.c_ijk);
	}
};


class FormatMatrixAsDiagonalAndComplement {
	exvector diagonal;
	string rest;
	static ex e_ij(int i,int j) {
		VectorField e_i("e_"+std::to_string(i+1)), e_j("e^"+std::to_string(j+1));
		return TensorProduct<VectorField,VectorField>(e_j,e_i);
	}	
public:
	FormatMatrixAsDiagonalAndComplement (const matrix& m) {		
		ex r;
		for (int i=0;i<m.rows();++i)
		for (int j=0;j<m.cols();++j)
			if (i==j) diagonal.push_back(m(i,i));
			else if (m(i,j)!=0) r+=collapse_roots(m(i,j))*e_ij(i,j);
		
		stringstream s;
		s<<latex;
		s<<to_latex_canonical_string(r);
		rest=s.str();
		if (rest[0]!='-') rest="+"+rest;	
	}
	string canonical_diagonal() const {
		vector<string> s;
		for (auto x: diagonal) s.push_back(to_latex_canonical_string(x));
		return horizontal(s);
	}
	string to_string() const {
		stringstream s;
		s<<"\\begin{array}{cc}"<<endl;		
		s<<"("<<canonical_diagonal()<<")\\\\"<<endl;
		s<<rest<<endl;
		s<<"\\end{array}"<<endl;
		return s.str();
	}
};



class PartitionedGeneralizedNilsolitonSurjectiveType {
	const vector<GeneralizedNilsolitonSurjectiveType> v;
	using Iterator = decltype(v.begin());
	map<HashGeneralizedNilsolitonSurjectiveTypeNoSignature,vector<Iterator>> partitions;
	HashGeneralizedNilsolitonSurjectiveTypeNoSignature hash_from_nilsoliton(const GeneralizedNilsolitonSurjectiveType& n) {
		auto hash=n.hash();
		return HashGeneralizedNilsolitonSurjectiveTypeNoSignature{hash.D, hash.c_ijk};
	}	
	void print_nilsolitons_and_signatures(ostream& os, const vector<Iterator>& I, const string& diagram_name) const {
		auto& g=*I.front();
		os<<"\\midrule"<<endl;
		os<<diagram_name<<"&";
		auto str_const=structure_constants_as_strings(g.lie_algebra());
		os<<Splitter{str_const,",",140}.to_string()<<"&";		
		os<<FormatMatrixAsDiagonalAndComplement{g.derivation()}.to_string()<<"&";
		vector<string> metrics;		
		for (auto i: I)
			metrics.push_back(horizontal(one_based_negative_indices(i->metric()),""));
		os<<Splitter{metrics,", ",40}.to_string_enclosed_in_curly_braces();
		os<<"\\\\"<<endl;
	}	

	void print_einstein_solvmanifold(ostream& os, const GeneralizedNilsolitonSurjectiveType& n) const {
		SemidirectProduct gtilde{n.lie_algebra(),n.derivation()};	
		os<<"--lie-algebra \""<<lie_group_to_string(gtilde)<<"\"";
		os<<" --indefinite-on-frame "<<horizontal(consecutive_numbers(1,gtilde.Dimension()+1));
		os<<" --timelike "<<horizontal(one_based_negative_indices(n.metric())," ")<<endl;	
	}
	void print_solvmanifolds_and_signatures(ostream& os, const vector<Iterator>& I) const {
		for (auto i: I)
			print_einstein_solvmanifold(os,*i);
	}
public:
	PartitionedGeneralizedNilsolitonSurjectiveType(const vector<GeneralizedNilsolitonSurjectiveType>& unpartitioned) : v{unpartitioned} {        
		for (auto i=v.begin();i!=v.end();++i)
			partitions[hash_from_nilsoliton(*i)].push_back(i);        
	}
	void print_nilsolitons(ostream& os, const string& diagram_name) const {
		for (auto& p: partitions) 
			print_nilsolitons_and_signatures(os,p.second,diagram_name);
	}	
	void print_einstein_solvmanifolds(ostream& os) const {
		for (auto& p: partitions) 
			print_solvmanifolds_and_signatures(os,p.second);
	}
};

