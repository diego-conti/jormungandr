exvector simplify_equation_with_sqrt(ex eqn) {
	eqn= eqn.subs(sqrt(wild(0))*sqrt(wild(1))==sqrt(wild(0)*wild(1)),subs_options::algebraic);	
	eqn=eqn.subs(sqrt(wild(0))-sqrt(wild(1))==abs(wild(0))-abs(wild(1)),subs_options::algebraic);	
	eqn=eqn.subs(abs(wild(0))*abs(wild(1))==abs(wild(0)*wild(1)),subs_options::algebraic);	
	eqn=eqn.expand();
	exmap matches;
	exvector result;	
	if (eqn.match(sqrt(wild(0))+sqrt(wild(1)),matches) ||eqn.match(-sqrt(wild(0))-sqrt(wild(1)),matches))
		result= {matches[wild(0)], matches[wild(1)]};
	else result={eqn};		
	for (auto& v: result) {
		exmap matches;		
		if (v.match(abs(wild(0)),matches))
			v=matches[wild(0)];								
	}
	return result;	
}

template<typename Iter1, typename Iter2>
exvector simplify_equations_with_sqrt(Iter1 begin, Iter2 end) {
	exvector result;
	for (auto i=begin;i!=end;++i) {
		auto simplified=simplify_equation_with_sqrt(*i);
		result.insert(result.end(),simplified.begin(),simplified.end());
	}
	return result;
}

ex collapse_roots(ex x) {
		ex subs=sqrt(wild(0))*sqrt(wild(1))==sqrt(wild(0)*wild(1));
		ex collapsed=x.subs(subs,subs_options::algebraic);
		while (x!=collapsed) {
			x=collapsed;
			collapsed=x.subs(subs,subs_options::algebraic);
		}
		return collapsed;
}
vector<string> structure_constants_as_strings(const LieGroup& G) {
		exvector s=G.StructureConstants();
		vector<string> r;		
		for (auto x: s) r.push_back(to_latex_canonical_string(collapse_roots(NormalizeRoots(x))));
		return r;
}

string lie_algebra_to_string(const LieGroup& G) {
		return horizontal(structure_constants_as_strings(G));
}