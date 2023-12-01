ex solve_wrt_symbols(const lst& eqns) {
	list<ex> unknowns;
	GetSymbols<symbol>(unknowns,eqns.begin(),eqns.end());	
	return lsolve(eqns,lst{unknowns});
}

//represents a derivation of the form a_1e^{i_1}\otimes e_{j_1}+... where the indices i_1,j_1,... are distinct
class NondiagonalMatrix {
	matrix d;
	map<int,int> i_to_j;
	static bool map_contains_pair(const map<int,int>& m, pair<int,int> p) {
		auto it=m.find(p.first);
		if (it==m.end()) return false;
		return it->second==p.second;
	}
	map<pair<int,int>,ex> symbols_;
	map<pair<int,int>,ex> create_symbols() const {
		map<pair<int,int>,ex> result;
		for (auto p: i_to_j) 
			result[p]=Parameter(N.A(p.second+1,p.first+1));
		return result;
	}
	NondiagonalMatrix(const matrix& m,const map<int,int>& nonzero) : d{m}, i_to_j{nonzero}, symbols_{create_symbols()} {
	}
public:
	static optional<NondiagonalMatrix> create_from_generic_element_and_nonzero_entries(const matrix& m, const map<int,int>& nonzero)  {
		lst eqns;
		for (int i=0;i<m.rows();++i)		
		for (int j=0;j<m.cols();++j) 						
			if (i!=j && !m(j,i).is_zero() && !map_contains_pair(nonzero,make_pair(i,j)))
				eqns.append(m(j,i).expand()==0);
		auto sol=solve_wrt_symbols(eqns);
		matrix part=ex_to<matrix>(ex(m).subs(sol));
		for (auto p: nonzero) 
			if (part(p.second,p.first).is_zero()) return nullopt;
		return NondiagonalMatrix{part,nonzero};
	}
	matrix as_matrix() const {return d;}
	auto pairs_begin() const {return i_to_j.begin();}
	auto pairs_end() const {return i_to_j.end();}
	map<pair<int,int>,ex> symbols() const {
		return symbols_;
	}
	bool operator<(const NondiagonalMatrix& other) const {return i_to_j<other.i_to_j;}
	string list_of_entries() const {
		stringstream s;		
		for (auto p: i_to_j) s<<p.first<<","<<p.second<<" ";		
		return s.str();
	}
};

class NondiagonalMatrixUptoAutomorphism {
	NondiagonalMatrix m;
	list<vector<int>> automorphisms_preserving_A;
	static map<pair<int,int>,ex> action(const vector<int>& automorphism, const map<pair<int,int>,ex>& symbols)  {
		map<pair<int,int>,ex> result;
		for (auto& p: symbols) {
			auto ij=p.first;
			auto f_ij=make_pair(automorphism[ij.first],automorphism[ij.second]);
			result[f_ij]=p.second;			
		}			
		return result;
	}
	static bool match(const map<pair<int,int>,ex>& s,const map<pair<int,int>,ex>& t) {
		if (s.size()!=t.size()) return false;
		return all_of(s.begin(),s.end(),[&t] (auto pair) {return t.find(pair.first)!=t.end();});
	}
	static vector<map<pair<int,int>,ex>> orbit(const map<pair<int,int>,ex>& symbols,  const list<vector<int>>& automorphisms) {
		vector<map<pair<int,int>,ex>> result{symbols};
		for (auto& x: automorphisms)
			result.push_back(action(x,symbols));
		return result;
	}
	static list<vector<int>> automorphisms_preserving(const map<pair<int,int>,ex>& s, const list<vector<int>>& automorphisms, int n) {
		list<vector<int>> result;
		vector<int> id(n);; 
		iota(id.begin(),id.end(),0);
		result.push_back(id);
		auto stabilizes_s = [&s] (auto x) {return match(s,action(x,s));};
		std::copy_if(automorphisms.begin(),automorphisms.end(),back_inserter(result),stabilizes_s);
		return result;
	}
public:
	NondiagonalMatrixUptoAutomorphism(const NondiagonalMatrix& m, const list<vector<int>>& automorphisms_preserving_A)
		: m{m}, automorphisms_preserving_A{automorphisms_preserving_A} {}
	NondiagonalMatrixUptoAutomorphism()=default;
	NondiagonalMatrixUptoAutomorphism(const NondiagonalMatrixUptoAutomorphism&)=default;
	operator NondiagonalMatrix() const {return m;}
	static vector<NondiagonalMatrixUptoAutomorphism> extract_section(list<NondiagonalMatrix> v, const list<vector<int>>& non_trivial_automorphisms) {
		constexpr bool DO_NOT_FACTOR=false;	//for debugging purposes
		vector<NondiagonalMatrixUptoAutomorphism> result;
		while (!v.empty()) {
			auto& i=v.front();
			result.emplace_back(i,automorphisms_preserving(i.symbols(),non_trivial_automorphisms,i.as_matrix().rows()));
            if (DO_NOT_FACTOR) v.pop_front();
			else {
	    		for (auto& x: orbit(i.symbols(),non_trivial_automorphisms)) {
		    		auto matches_x=[&x] (auto d) {return match(x,d.symbols());};
			    	v.remove_if(matches_x);				
			    }
            }
		}
		return result;
	}
	list<vector<int>> stabilizer() const {return automorphisms_preserving_A;}
};
