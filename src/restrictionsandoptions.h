//class defining restrictions on the sets of indices A that should be considered
struct Restrictions {		
	EqualityOrInequality condition_on_size_of_A,condition_on_free_parameters;

	int no_free_parameters(const map<int,int>& A, const GL& gl, const Derivations& der)  const {
		auto diagonal_derivations=der.diagonal_derivations();
		auto gen_diag=gl.glToMatrix(diagonal_derivations.GenericElement());
		lst eqns;
		for (auto& aij: A)	{
			int i=aij.first, j=aij.second;
			eqns.append(gen_diag(i,i)-gen_diag(j,j)-gen_diag.trace());
		}				
		auto subspace=diagonal_derivations.SubspaceFromEquations(eqns.begin(),eqns.end());
		auto codimension=diagonal_derivations.Dimension()-subspace.Dimension();
		return A.size()-codimension;
	}
public:
	Restrictions free_parameters(string condition) const {
		auto result=*this;
		result.condition_on_free_parameters=EqualityOrInequality{condition};
		return result;
	}
	Restrictions bound_on_A(string bound_on_A) const {
		auto result=*this;
		result.condition_on_size_of_A=EqualityOrInequality{bound_on_A};
		return result;
	}
	bool satisfied(const map<int,int>& A, const GL& gl, const Derivations& der) const {
		return condition_on_size_of_A.verified_by(A.size()) &&		
			condition_on_free_parameters.verified_by(no_free_parameters(A,gl,der));
	}
};

struct OutputOptions {
	shared_ptr<ostream> metrics_stream = make_shared<ostream>(nullptr);
	bool verbose=false;
};