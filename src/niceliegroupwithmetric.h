//a class for nice lie groups where the structure constants are determined by the vector X in such a way that the Ricci tensor is forced to take a specific form

class CoefficientConfigurationWithAssignedX : public CoefficientConfigurationWithVariableSigns {
	static Weights weights_and_values(const WeightBasis& weight_basis, const exvector& X) {
		Weights weights;
		auto coeff=X.begin();
		nice_log<<weight_basis.weights_and_coefficients()<<endl;
		for (auto weight: weight_basis.weights_and_coefficients()) 
			weights.add_weight(weight,sqrt(abs(*coeff++)));  
		return weights;	
	}

public:
	CoefficientConfigurationWithAssignedX(const WeightBasisAndProperties& weight_basis, const exvector& X) 
		: CoefficientConfigurationWithVariableSigns{weight_basis.sign_configurations(),weight_basis.number_of_nodes(),weights_and_values(weight_basis,X)} {}
};


class NiceLieGroupWithMetric : public LieGroupsFromDiagram {
	NiceLieGroupWithMetric(const CoefficientConfigurationWithAssignedX& configuration)	
		: LieGroupsFromDiagram(configuration.lie_algebra_dimension()) {	 
		ExVector de(configuration.lie_algebra_dimension());
		int i=0;
		for (WeightAndValue weight : configuration.weights()) 
			de[weight.node_out]+=weight.value*e()[weight.node_in1]*e()[weight.node_in2];
		for (int i=1;i<=de.size();++i)
			Declare_d(e(i), de(i));
	}
	static	list<NiceLieGroupWithMetric> from_coefficient_configuration (CoefficientConfigurationWithAssignedX&& configuration) {
		list<NiceLieGroupWithMetric> result;
		insert_new_lie_group(result,configuration);
		while (configuration) {
			insert_new_lie_group(result,configuration);
			++configuration;
		}
		return result;
	}
	static void insert_new_lie_group(list<NiceLieGroupWithMetric>& out_list, const CoefficientConfigurationWithAssignedX& configuration) {
		NiceLieGroupWithMetric nice_lie_group{configuration};
		if (nice_lie_group.solve_linear_ddzero()) {		
			for (auto x: nice_lie_group.e()) {
				auto ddx=nice_lie_group.d(nice_lie_group.d(x));
				if (!ddx.is_zero()) nice_log<<"d^2!=0 in general, "<<ddx<<endl;
			}
			out_list.push_back(move(nice_lie_group));
		}
	}
public:
	static list<NiceLieGroupWithMetric> from_weight_basis_and_X(const WeightBasisAndProperties& weight_basis, const exvector& X) {
		CoefficientConfigurationWithAssignedX configuration{weight_basis,X};
    	return from_coefficient_configuration(move(configuration));    	
	}
};