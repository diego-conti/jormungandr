//return a semidirect product lie algebra, imposing the following:
//d\tilde e^i = de^i -D_ij e^0\wedge e^j   
//e^i[e_0,e_j]=e^i(De_j) = D(i,j)
class SemidirectProduct : public LieGroupHasParameters<true>, public ConcreteManifold {
public:
	SemidirectProduct(const LieGroupHasParameters<true>& G, const matrix& D) : ConcreteManifold(G.Dimension()+1) {
		lst subs=LinearMapToSubstitutions(G.e(),e().begin(),e().begin()+G.Dimension());	
		ex e0=e(Dimension());
		Has_dTable::Declare_d(e0,0);
		for (int i=0;i<G.Dimension();++i) {
			ex D_i;
			for (int j=0;j<G.Dimension();++j)
				D_i+=D(i,j)*e()[j];			
			Has_dTable::Declare_d(e()[i], NormalizeRoots(G.d(G.e()[i]).subs(subs)-e0*D_i));
		}

	}
};
