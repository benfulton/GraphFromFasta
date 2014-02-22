#include <string>

class NonRedKmerTable;
class DNAVector;
class config;
class vecDNAVector;

class Welder
{
public:
    Welder(int k, int kk);
    
    void SetTable(NonRedKmerTable * p);

    // constructs the required weldable kmer from two contigs a,b and positions one and two
    void WeldableKmer(DNAVector & out, 
                      const DNAVector & a, int one, 
                      const DNAVector & b, int two); 
    
    bool Weldable(const DNAVector & a, int one, const DNAVector & b, int two, int thresh, std::string& welding_kmer, int& weldable_kmer_read_count) ;

    int weldmer_count_in_reads(const DNAVector weldmer);

    
private:
    NonRedKmerTable * m_pTab;
    
    int m_k;
    int m_kk;
};

bool SimpleHalves(const DNAVector & d, float MIN_WELD_ENTROPY, float MAX_RATIO_INTERNALLY_REPETITIVE);
void Add(vecDNAVector & all, DNAVector & add, int & counter) ;
	