//write a sequence of latex strings possibly separating into multiple lines
class Splitter {
	const vector<string> seq;
	const string sep;
	const int max_size;
	vector<string> lines;
	vector<string>::const_iterator fill_line(vector<string>::const_iterator i) {
		stringstream ss;
		ss<<*i;
		int size=0;
		while (++i!=seq.end()) { 
			size+=i->size()+sep.size();
			if (ss.str().size() && size>max_size) break;
			ss<<sep<<*i;
		}
		lines.push_back(ss.str());
		return i;
	}
	void print_lines_in_array(ostream& os) const {
		auto i=lines.begin();
		os<<*i;
		while (++i!=lines.end()) 
			os<<"\\\\"<<endl<<*i;
	}
public:
	Splitter(const vector<string>& sequence, string sep, int max_size=15) : seq{sequence}, sep{sep},max_size{max_size} {
		auto i=seq.begin();
		while (i!=seq.end()) {
			if (!lines.empty()) lines.back()+=sep;
			i=fill_line(i);
		}
	}
	string to_string(string open="", string close="") const {
		if (lines.size()==0) return {};
		else if (lines.size()==1) return open+lines[0]+close;
		stringstream ss;
		ss<<"\\begin{array}{c}"<<endl;
		ss<<open;		
		print_lines_in_array(ss);
		ss<<close<<endl;
		ss<<"\\end{array}"<<endl;
		return ss.str();
	}
	string to_string_enclosed_in_curly_braces() const {
		return to_string("\\{","\\}");
	}
};