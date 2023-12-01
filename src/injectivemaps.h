
class PartialMapIterator {
    map<int,int> m;
    int n;

    static bool pairs_overlap(pair<int,int> a, pair<int,int> b) {
        return a.first==b.first || a.second==b.first || a.first==b.second || a.second==b.second;
    }
    bool injective_on_complement_up_to(map<int,int>::iterator i) const {
        if (i->first==i->second) return false;
        for (auto j=m.begin();j!=i;++j)
            if (pairs_overlap(*i,*j)) return false;
        return true;
    }
    //initialize elmeents after i in such a way that they give the first admissible map in lexicographic order
    bool set_first_map_after(map<int,int>::iterator i) {                
        while (++i!=m.end()) {                        
            m[i->first]=0; 
            while (m[i->first]<n && !injective_on_complement_up_to(i)) 
                ++m[i->first];
            if (m[i->first]==n) return false;            
        }        
        return true;
    }
    optional<int> first_element_not_in_domain() const {
        for (int i=0;i<n;++i) 
            if (m.find(i)==m.end()) return i;
        return nullopt;
    }
    //choose first map with the given domain
    bool set_first_map() {
        if (n<1) return false;
        int k=m.begin()->first;
        auto first=first_element_not_in_domain();
        if (first) {
            m[k]=first.value();
            return set_first_map_after(m.begin());
        }
        else return false;
    }

    bool increment_at(map<int,int>::iterator i) {
        while (++m[i->first]!=n) 
            if (injective_on_complement_up_to(i) && set_first_map_after(i)) return true;
        return false;    
    }
    bool increment_map() {
        auto i=m.end();
        while (i!=m.begin()) 
            if (increment_at(--i)) return true;
        return false;
     }
     //iterate over subsets like this: choose the maximum element k in the set; if k=n-1, remove it and increment the maximum; otherwise add k+1 to the set
    bool increment_domain() {
        do {
            int k=m.rbegin()->first;
            if (k==n-1) {
                m.erase(k);
                if (m.empty()) return false;
                k=m.rbegin()->first;
                m.erase(k);
            }
            m[k+1]=0;   
        }
        while (!set_first_map());
        return true;
    }

    PartialMapIterator(int n) : n{n} {}
public:
    static PartialMapIterator begin(int n) {    
        PartialMapIterator i{n};
        i.m[0]=1;
        return i;
    }
    static PartialMapIterator end(int n) {return PartialMapIterator{n};}
    bool operator!=(const PartialMapIterator& i) {return m!=i.m;}
    void operator++() {        
        if (!increment_map() && !increment_domain()) m.clear();
    }
    map<int,int> operator*() const {return m;}
};

//given a finite set {x_1,..,x_n}, iterate through partial injective maps that map some subset {x_1,..,x_k} to the complement {x_k+1,...,x_n}
class PartialInjectiveMapsToComplement {
    int n;
public:
    PartialInjectiveMapsToComplement(int n) : n{n} {}
    PartialMapIterator begin() const {return PartialMapIterator::begin(n);}
    PartialMapIterator end() const {return PartialMapIterator::end(n);}
};
