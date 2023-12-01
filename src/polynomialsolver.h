class Solution {	
	ex solution;
public:
	Solution(ex solution) : solution{solution} {}
	ex substitute(ex x) const {
		return x.subs(solution);
	}
	exvector substitute(exvector x) const {
		transform(x.begin(),x.end(),x.begin(),[this] (ex x) {return x.subs(solution);});
		return x;
	}
	template<typename T>
    map<T,ex> substitute(map<T,ex> map_to_ex) const {
	    for (auto i= map_to_ex.begin();i!=map_to_ex.end();++i)
            i->second=i->second.subs(solution);
		return map_to_ex;
	}
	matrix substitute(matrix x) const {
		for (int i=0;i<x.nops();++i)
			x.let_op(i)=x.op(i).subs(solution);
		return x;
	}
	template<typename T, typename U=typename std::enable_if_t<std::is_base_of<LieGroup, T>::value>>
	T substitute(T lie_algebra) {		
		lst as_lst= (is_a<lst>(solution))? ex_to<lst>(solution) : lst{solution};
		lie_algebra.DeclareConditions(as_lst);
		return lie_algebra;
	}
};

class SolutionsIterator {
	exvector::const_iterator iter;
public:
	SolutionsIterator(exvector::const_iterator i) : iter{i} {}
	SolutionsIterator operator++() {++iter; return *this;}
	bool operator!=(SolutionsIterator other) {return iter!=other.iter;}
	Solution operator*() const {return *iter;}	
};

class Solutions {
	exvector solutions_;
public:
	template<typename Iter1, typename Iter2> Solutions(Iter1 begin, Iter2 end) : solutions_{begin,end} {}
	SolutionsIterator begin() const {return solutions_.begin();}
	SolutionsIterator end() const {return solutions_.end();}
};

template<typename Unknown> 
class PolynomialSolver {
	exvector solutions_;	
	bool can_solve_;
	void solve_no_variables(ex equation) {
		can_solve_=true;
		if (equation.is_zero())
			solutions_.push_back(lst{});
	}

    void solve_in_variable(ex equation, ex x) {
        exmap matches;
        if (equation.match(abs(wild(0))-abs(wild(1)),matches)) {
            solve_polynomial(matches[wild(0)]-matches[wild(1)],x);
            solve_polynomial(matches[wild(0)]+matches[wild(1)],x);
        }
        else solve_polynomial(equation,x);
    }
	void solve_polynomial(ex equation, ex x) {    
        if (!equation.is_polynomial(x)) {
            can_solve_=false; return;
        }
        equation=equation.expand();    
		auto a=equation.coeff(x,2);
		auto b=equation.coeff(x,1);
		auto c=equation.coeff(x,0);
		auto delta=b*b-4*a*c;
		switch (equation.degree(x)) {
			case 1: 
				can_solve_=true;
				solutions_.push_back(x==-c/b);
				break;
			case 2: 
                assert(a!=0);
				can_solve_=true;
				if (delta.is_zero()) solutions_.push_back(x==-b/(2*a));
				else if (delta>0) {
					solutions_.push_back(x==(-b+sqrt(delta))/(2*a));
					solutions_.push_back(x==(-b-sqrt(delta))/(2*a));
				}
				break;
			case 0:				
				assert(false);
				break;
			default:
				can_solve_=false;
			}
	}	
	PolynomialSolver()=default;
public:	
	PolynomialSolver(ex equation) {
		list<ex> parameters;		
		GetSymbols<StructureConstant> (parameters,equation);
		if (parameters.empty()) {
			list<ex> variables;		
			GetSymbols<Unknown> (variables,equation);
			if (variables.empty()) solve_no_variables(equation);
			else if (variables.size()>1) can_solve_=false;
			else solve_in_variable(equation,variables.front());
		}
		else can_solve_=false;
	}
	static PolynomialSolver solver_that_cannot_solve() {
		PolynomialSolver result;
		result.can_solve_=false;
		return result;	
	}
	bool can_solve() const {return can_solve_;}
	bool has_solution() const {return solutions_.size();}
	list<ex> substitute_solutions(ex expression) const {
		list<ex> result;
		transform(solutions_.begin(),solutions_.end(),back_inserter(result),[expression] (ex solution) {return expression.subs(solution);});
		return result;
	}
	list<exvector> substitute_solutions(const exvector& expression) const {
		auto subs_into_vector=[&expression] (ex solution) {
			exvector result;
			transform(expression.begin(),expression.end(),back_inserter(result),[solution] (ex x) {return x.subs(solution);});
			return result;
		};	
		list<exvector> result;
		transform(solutions_.begin(),solutions_.end(),back_inserter(result),subs_into_vector);
		return result;
	}

    template<typename T>
    list<map<T,ex>> substitute_solutions(const map<T,ex>& map_to_ex) const {
		auto subs_into_map=[&map_to_ex] (ex solution) {
			map<T,ex> result;
            for (auto p : map_to_ex)
                result[p.first]=p.second.subs(solution);
			return result;
		};	
		list<map<T,ex>> result;
		transform(solutions_.begin(),solutions_.end(),back_inserter(result),subs_into_map);
		return result;

    }
	void add_solutions_from(const PolynomialSolver& solver) {
		if (!solver.can_solve()) can_solve_=false;
		solutions_.insert(solutions_.end(),solver.solutions_.begin(),solver.solutions_.end());
	}	

	Solutions solutions() const {
		return Solutions{solutions_.begin(),solutions_.end()};
	}
};

