
exvector X_solving_generalized_nilsoliton(const WeightMatrix& weight_matrix, const ExVector& im_X) {		
	auto MDelta=weight_matrix.M_Delta();	
	auto complete_matrix = adjoin(transpose(MDelta),ColumnVectorMatrixBuilder{im_X});
	auto X=	solve_over_Q(complete_matrix,generate_variables<Unknown>(N.x,weight_matrix.rows()));
	return X;
};


exvector X_for_generalized_nilsolitons(const LabeledTree& diagram, const WeightMatrix& M_Delta, map<pair<int,int>,ex> A, const matrix& D, const NondiagonalMatrix& d) {	
	ex sum=std::transform_reduce(A.begin(),A.end(),ex{0}, std::plus{}, [](auto p) {return p.second;});
	ExVector im_X(diagram.number_of_nodes(),-1-sum);
	for (auto p : A) {
		auto ij =p.first;
		auto A_ji=p.second;
		im_X[ij.first]+=A_ji;
		im_X[ij.second]-=A_ji;		
	}	
	for (int i=0;i<im_X.size();++i)
		im_X[i]+=D(i,i);	
	for (auto& x: im_X) x*=2*D.trace();
	return X_solving_generalized_nilsoliton(M_Delta,im_X);
}

//return the equation Tr (DX)==(1+\sum A_ij)Tr X+\sum A_{ij}(mu_j-mu_i)
ex get_generalized_nik_condition(ex D, ex X, const NondiagonalMatrix& d) {
	auto A=d.symbols();
	auto tr_XD=ex_to<matrix>((X*D).evalm()).trace();
	auto tr_X=ex_to<matrix>(X).trace();
	ex RHS=tr_X;
	for (auto it=d.pairs_begin();it!=d.pairs_end();++it) {
		auto i=it->first;
		auto j=it->second;
		auto mu_j_minus_mu_i=ex_to<matrix>(X)(j,j)-ex_to<matrix>(X)(i,i);
		RHS+=A.at(make_pair(i,j))*(tr_X+mu_j_minus_mu_i);
	}
	return tr_XD==RHS.expand();
}

//if the diagonal part of D is (lambda_1,...,lambda_n), return the equations lambda_1+...+\lambda_n == (lambda_i-lambda_j), where (i,j) ranges in A
lst trace_D_equals_weights (const matrix& D, const NondiagonalMatrix& d) {
	lst result;
	for (auto it=d.pairs_begin();it!=d.pairs_end();++it) {
		auto i=it->first;
		auto j=it->second;
		result.append(D.trace()==D(i,i)-D(j,j));
	}
	return result;
}


//return zero if e^j\otimes e_i\circ \ad _v has trace zero for all v, which means that e_i\hook de^j=0
bool tr_ad_v_circ_Dstar_is_zero(const NiceLieGroup& G, const NondiagonalMatrix& d) {	
	auto tr_ad_v_circ_Dstar_is_zero = [&G] (auto pair) {return 	Hook(G.e()[pair.first], G.d(G.e()[pair.second])).is_zero();};
	return all_of(d.pairs_begin(),d.pairs_end(),tr_ad_v_circ_Dstar_is_zero);
}

class NondiagonalGeneralizedNilsolitons {
	const LabeledTree& diagram;
	WeightMatrix weight_matrix;
	GL gl;
	VectorSpace<DifferentialForm> diagonal_derivations;
	list<NondiagonalMatrixUptoAutomorphism> all_nondiagonal_derivations;
	ostream& stream_with_metrics_found;
	OutputOptions options;
	list<NondiagonalMatrix> nondiagonal_derivations(const NiceLieGroup& G, const Derivations& der, Restrictions restrictions) const {				
		auto offdiag=der.space_containing_offdiagonal_derivations();
		auto D=gl.glToMatrix(offdiag.GenericElement()+diagonal_derivations.GenericElement());
		list<NondiagonalMatrix> result;
		for (auto x: PartialInjectiveMapsToComplement(G.Dimension())) {		
			if (!restrictions.satisfied(x,gl,der)) continue;			
			auto der=NondiagonalMatrix::create_from_generic_element_and_nonzero_entries(D,x);			
			if (der && tr_ad_v_circ_Dstar_is_zero(G,der.value())) {			
				result.push_back(der.value());
			}
		}
		return result;
	}
	void impose_generalized_nik_condition(lst& eqns, ex D, const NondiagonalMatrix& d) const {	
		for (auto f: diagonal_derivations.e()) 
			eqns.append(::get_generalized_nik_condition(D,gl.glToMatrix(f),d));	
	}
	static list<matrix> assign_offdiagonal_entries(matrix D, pair<int,int> ij, ex A_ij, const map<pair<int,int>,ex>& rest) {
		ex tr_D=D.trace();
		D(ij.second,ij.first)=NormalizeRoots(sqrt(abs(A_ij*2*tr_D)));
		auto with_plus=assign_offdiagonal_entries(D,rest);
		D(ij.second,ij.first)=-D(ij.second,ij.first);
		auto with_minus=assign_offdiagonal_entries(D,rest);		
		with_plus.splice(with_plus.end(),with_minus);
		return with_plus;
	}

	static list<matrix> assign_offdiagonal_entries(const matrix& D, const map<pair<int,int>,ex>& A) {
		if (A.empty()) return {D};		
		auto a=*A.begin();	
		auto rest=A;
		rest.erase(rest.begin());
		auto ij=a.first;
		auto A_ij=a.second;
		return assign_offdiagonal_entries(D,ij,A_ij, rest);
	}
public:
	NondiagonalGeneralizedNilsolitons(const LabeledTree& diagram,  Restrictions restrictions={}, OutputOptions options={}) :
		 diagram{diagram}, weight_matrix{diagram.weights(),diagram.number_of_nodes()}, gl{diagram.number_of_nodes()},stream_with_metrics_found{*options.metrics_stream}, options{options} {
		auto lie_algebras=NiceLieGroup::from_weight_basis(diagram.weight_basis(processor())); 	
		for (auto lie_algebra : lie_algebras) {
			Derivations der{gl,lie_algebra};
			diagonal_derivations=der.diagonal_derivations();
			auto nondiagonal=nondiagonal_derivations(lie_algebra,der,restrictions);	//compute elements \mathcal A satisfying N1-N3			
			auto section=NondiagonalMatrixUptoAutomorphism::extract_section(nondiagonal,diagram.nontrivial_automorphisms());
			all_nondiagonal_derivations.insert(all_nondiagonal_derivations.end(),section.begin(),section.end());			
		}
	}
	vector<GeneralizedNilsoliton> generalized_nilsolitons_from_nondiagonal_derivation(const NondiagonalMatrix& derivation) const {		
		vector<GeneralizedNilsoliton> result;
		auto D=derivation.as_matrix();
		lst eqns=trace_D_equals_weights(D,derivation);
		auto A=derivation.symbols();
		impose_generalized_nik_condition(eqns,D,derivation);
		//impose N4 and N5
		auto sol=solve_wrt_symbols(eqns);					
		auto derivation_compatible_with_generalized_nik=ex_to<matrix>(ex(D).subs(sol));				
		
		for (auto& a: A) {
			a.second=a.second.subs(sol);			
		}
		
		//compute X so that equation (4) holds
		auto X=X_for_generalized_nilsolitons(diagram,weight_matrix,A,derivation_compatible_with_generalized_nik,derivation);			
		for (auto derivation : assign_offdiagonal_entries(derivation_compatible_with_generalized_nik,A)) 
			result.emplace_back(derivation,A,X);
		return result;
	}

	//return the generalized nilsolitons satisfying the restrictions if the root matrix is surjective, or the empty set otherwise.
	vector<GeneralizedNilsolitonSurjectiveType> generalized_nilsolitons_surjective_type() const {
		if (weight_matrix.rank_over_Q()!=weight_matrix.rows()) return {};
		vector<GeneralizedNilsolitonSurjectiveType> result;
		for (auto& x: generalized_nilsolitons()) {
			auto n=GeneralizedNilsolitonSurjectiveType::expand(x,diagram,weight_matrix,gl,options);
			result.insert(result.end(),n.begin(),n.end());
		}
		return result;
	}

	//return a section for the generalized nilsolitons satisfying the restrictions if the root matrix is surjective, or the empty set otherwise.
	vector<GeneralizedNilsolitonSurjectiveType> section_of_generalized_nilsolitons_surjective_type() const {
		if (weight_matrix.rank_over_Q()!=weight_matrix.rows()) return {};
		vector<GeneralizedNilsolitonSurjectiveType> result;
		for (auto& x: all_nondiagonal_derivations) {
			auto v=generalized_nilsolitons_from_nondiagonal_derivation(x);
			auto n=GeneralizedNilsolitonSurjectiveType::expand(v,diagram,weight_matrix,gl,options);			
			auto f=factor_out_by(n, x.stabilizer());
			result.insert(result.end(),f.begin(),f.end());
		}
		return result;
	}

	//return the generalized nilsolitons satisfying the restrictions
	vector<GeneralizedNilsoliton> generalized_nilsolitons() const {
		vector<GeneralizedNilsoliton> result;
		for (auto x: all_nondiagonal_derivations) {
			auto n=generalized_nilsolitons_from_nondiagonal_derivation(x);			
			result.insert(result.end(),n.begin(),n.end());
		}
		return result;
	}
};
