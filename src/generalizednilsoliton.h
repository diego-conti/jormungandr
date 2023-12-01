#include "semidirectproduct.h"
#include "sqrt.h"

DiagramProcessor processor() {
	DiagramProcessor processor{with_lie_algebra};
	Filter filter;
	//filter.only_cokernel_MDelta_dimension("=0"s);
	processor.setFilter(filter);	
	processor.set(DiagramDataOption::with_diagonal_nilsoliton_metrics);	
	processor.invert_nodes();	
	return processor;	
}


WEDGE_DECLARE_NAMED_ALGEBRAIC(Parameter,symbol)

struct GeneralizedNilsoliton {
	matrix D;
	map<pair<int,int>,ex> A;
	exvector X;
	GeneralizedNilsoliton(const matrix& D, const map<pair<int,int>,ex>& A, const exvector& X) : D{D},A{A},X{X} {}
	void print(ostream& os) const {
		os<<"D="<<D<<endl;
		os<<"X="<<horizontal(X)<<endl;
		for (auto p: A) {
			auto ij=p.first;
			auto A_ij=p.second;
			os<<"A_"<<ij.second+1<<ij.first+1<<"="<<A_ij<<endl;
		}
	}
};

namespace ActionOfAutomorphisms {
	 map<tuple<int,int,int>,ex> action(const vector<int>& automorphism , const SignConfiguration& delta,const map<tuple<int,int,int>,ex>& c_ijk)  {
		map<tuple<int,int,int>,ex> result;
		for (auto p: c_ijk) {
			int i=get<0>(p.first), j=get<1>(p.first), k=get<2>(p.first);
			auto f_i=automorphism[i], f_j=automorphism[j], f_k=automorphism[k];
			if (f_i>f_j) {swap(f_i,f_j); p.second=-p.second;}
			auto f_ijk=make_tuple(f_i,f_j,f_k);		 	
			result[f_ijk]=p.second*delta[i]*delta[j]*delta[k];			
		}			
		return result;
	}
	matrix action(const vector<int>& f, const SignConfiguration& delta, const matrix& m) {
		matrix fm(m.rows(),m.cols());
		for (int i=0;i<m.rows();++i)
		for (int j=0;j<m.cols();++j)
			fm(f[i],f[j])=m(i,j)*delta[i]*delta[j];
		return fm;
	}
	vector<int> action(const vector<int>& f, const SignConfiguration& delta, const vector<int>& signature) {
		vector<int> f_signature(signature.size());
		for (int i=0;i<signature.size();++i)			
			f_signature[f[i]]=signature[i];
		return f_signature;
	}
}


struct HashGeneralizedNilsolitonSurjectiveType {
	matrix D;	
	vector<int> signature;	
	map<tuple<int,int,int>,ex> c_ijk;	//the structure constants, i.e. [e_i,e_j]=c_ij e_k, where i<j 	

	HashGeneralizedNilsolitonSurjectiveType action(const vector<int>& f, const SignConfiguration& delta) const {
		using ActionOfAutomorphisms::action;
		return HashGeneralizedNilsolitonSurjectiveType{action(f,delta,D),  action(f,delta,signature),action(f,delta,c_ijk)};		
	}	
	bool operator==(const HashGeneralizedNilsolitonSurjectiveType& o) const {
		return D==o.D && signature==o.signature && c_ijk==o.c_ijk;
	}
	void print(ostream& os) const {
		os<<"D="<<D<<endl;
		os<<"signature="<<horizontal(signature)<<endl;
		for (auto p: c_ijk) {
			int i=get<0>(p.first)+1, j=get<1>(p.first)+1, k=get<2>(p.first)+1;
			os<<"[e_"<<i<<","<<"e_"<<j<<"]="<<p.second<<" e_"<<k<<endl;
		}
	}
};

class GeneralizedNilsolitonSurjectiveType {
	matrix D;
	map<pair<int,int>,ex> A;
	NiceLieGroupWithMetric g;
	vector<int> signature;	
	list<Weight> weights;
public:
	GeneralizedNilsolitonSurjectiveType(const matrix& D, const	map<pair<int,int>,ex>& A, const NiceLieGroupWithMetric& g, const vector<int>& signature, const list<Weight>& weights)
		: D{D},A{A},g{g},signature{signature},weights{weights} {			
		}
	static vector<GeneralizedNilsolitonSurjectiveType> expand(const vector<GeneralizedNilsoliton>& nil, const LabeledTree& diagram, const WeightMatrix& weight_matrix, const GL& gl, OutputOptions options) {
		 vector<GeneralizedNilsolitonSurjectiveType> result;
		 for (auto x: nil) {
			auto v=expand(x,diagram,weight_matrix,gl,options);
		 	result.insert(result.end(),v.begin(),v.end());
		 }
		 return result;
	}
	static vector<GeneralizedNilsolitonSurjectiveType> expand(const GeneralizedNilsoliton& nil, const LabeledTree& diagram, const WeightMatrix& weight_matrix, const GL& gl, OutputOptions options) {
		vector<GeneralizedNilsolitonSurjectiveType> result;
		auto& weight_basis=diagram.weight_basis(processor());
		for (auto& g: NiceLieGroupWithMetric::from_weight_basis_and_X(weight_basis,nil.X)) {
			lst derivation_if=equations_such_that_linear_map_is_derivation(g,gl,gl.MatrixTo_gl(nil.D));			
			auto simplified=simplify_equations_with_sqrt(derivation_if.begin(),derivation_if.end());
			if (options.verbose && !simplified.empty()) cout<<"equations: "<<horizontal(simplified)<<endl;
			solve_eqns_and_add(g,nil.A,nil.D,nil.X,weight_basis,simplified,weight_matrix,diagram.weights(),result);
		}
		return result;
	}
	void print(ostream& os) const {
		os<<"generalized nilsoliton of surjective type:"<<endl;
		os<<"\\g="<<horizontal(g.StructureConstants())<<endl;
		os<<"g=\\diag("<<horizontal(signature)<<")"<<endl;
		os<<horizontal(one_based_negative_indices(signature),",")<<endl;
		os<<"D="<<D<<endl;
		for (auto p: A) {
			auto ij=p.first;
			auto A_ij=p.second;
			os<<"A_"<<ij.second+1<<ij.first+1<<"="<<A_ij<<endl;
		}
	}
	HashGeneralizedNilsolitonSurjectiveType hash() const {
		return {D,signature,structure_constants(g)};
	}
	int Dimension() const {return g.Dimension();}
	matrix derivation() const {return D;}
	const NiceLieGroupWithMetric lie_algebra() const {return g;}
	SignConfiguration metric() const {return signature;}
private:
	 map<tuple<int,int,int>,ex> structure_constants(const LieGroupsFromDiagram& g) const {
		map<tuple<int,int,int>,ex> c_ijk;	//the structure constants, i.e. [e_i,e_j]=c_ij e_k, where i<j and e_k is determined by the nice condition.	
		auto c=g.c(weights);
		auto i=c.begin();
		for (auto w: weights) {
			assert(w.node_in1<w.node_in2);
			c_ijk[make_tuple(w.node_in1,w.node_in2,w.node_out)]=*i++;
		}
		return c_ijk;
	}
	static bool signs_compatible(const map<pair<int,int>,ex>& A, const matrix& D, const vector<int>&  signature)  {
		ex tr_D=D.trace();
		auto is_compatible = [tr_D,&signature] (auto p) {		
			int i=p.first.first;
			int j=p.first.second;
			auto A_ij=p.second;					
			return A_ij*signature[i]*signature[j]*tr_D>0;
		};
		return all_of(A.begin(),A.end(),is_compatible);
	}

	static void add_case_with_no_parameters(const NiceLieGroupWithMetric& g, const map<pair<int,int>,ex>& A, const matrix& D, const exvector& X, const WeightBasis& weight_basis,const WeightMatrix& weight_matrix,const list<Weight>& weights, vector<GeneralizedNilsolitonSurjectiveType>& v)  {					
			DiagonalMetric metric("generalized nilsoliton",weight_matrix,X);		
			for (auto& signature : metric.signatures(g.csquared(weight_basis)).as_vectors())
				if (signs_compatible(A,D,signature)) 
					v.emplace_back(D,A,g,signature,weights);
	}		

	static void solve_eqns_and_add(const NiceLieGroupWithMetric& g, const map<pair<int,int>,ex>& A, const matrix& D, const exvector& X, const WeightBasis& weight_basis, exvector equations, const WeightMatrix& weight_matrix,const list<Weight> weights,vector<GeneralizedNilsolitonSurjectiveType>& v)  {		
		if (!equations.size()) add_case_with_no_parameters(g,A,D,X,weight_basis,weight_matrix,weights,v);		
		else {			
			auto first=equations.back();
			auto rest=equations; 
			rest.pop_back();			
			PolynomialSolver<Parameter> solver(first);
			if (solver.can_solve()) 
				for (auto sol : solver.solutions()) {
					solve_eqns_and_add(sol.substitute(g), sol.substitute(A), sol.substitute(D),sol.substitute(X),weight_basis, sol.substitute(rest),weight_matrix,weights,v);
				}
			else {
				cout<<"cannot solve the equations automatically:"<<horizontal(equations,"; ")<<endl;
				cout<<"dflt form:"<<dflt<<equations<<latex<<endl;
				cout<<"on Lie algebra "<<horizontal(g.StructureConstants())<<endl;
				cout<<"with D="<<D<<endl;
				cout<<"X="<<horizontal(X)<<endl<<endl;
			}
		}
	}
};

template<typename T,typename F,typename G>
vector<T> extract_section(list<T> X, F orbit_of, G matches) {
	vector<T> result;
		while (!X.empty()) {
			auto& i=X.front();
			result.emplace_back(i);
			for (auto& x: orbit_of(i)) {
				auto matches_x=matches(x);
				X.remove_if(matches_x);
			}
		}
		return result;
}

	
vector<GeneralizedNilsolitonSurjectiveType> factor_out_by(list<GeneralizedNilsolitonSurjectiveType> v, const list<vector<int>>& automorphisms) {
	auto orbit_of=[&automorphisms] (const GeneralizedNilsolitonSurjectiveType& g) {
		HashGeneralizedNilsolitonSurjectiveType hash=g.hash();
		vector<HashGeneralizedNilsolitonSurjectiveType> hashes;
		for (auto& a: automorphisms) {
			SignConfiguration delta(g.Dimension());
			hashes.push_back(hash.action(a,delta));
			while (delta.has_next()) 
				hashes.push_back(hash.action(a,++delta));					
		} 		
		return hashes;
	};
	auto matches=[] (const HashGeneralizedNilsolitonSurjectiveType& hash) {		
		return [&hash] (const GeneralizedNilsolitonSurjectiveType& y) {
			return hash==y.hash();
		};
	};
	auto section= extract_section(std::move(v),orbit_of,matches);
	return section;
}

vector<GeneralizedNilsolitonSurjectiveType> factor_out_by(const vector<GeneralizedNilsolitonSurjectiveType>& v, const list<vector<int>>& automorphisms) {
	return factor_out_by(list<GeneralizedNilsolitonSurjectiveType>{v.begin(),v.end()},automorphisms);
}
