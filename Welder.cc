#include "Welder.h"
#include "DNAVector.h"
#include "NonRedKmerTable.h"
#include "DNAUtil.h"
#include "config.h"

#undef DEBUG

static bool DEBUG = false;


Welder::Welder(int k, int kk) {
    m_k = k;
    m_kk = kk;
    m_pTab = NULL;
}
    
void Welder::SetTable(NonRedKmerTable * p) {
    m_pTab = p;
}
    

// constructs the required weldable kmer from two contigs a,b and positions one and two
void Welder::WeldableKmer(DNAVector & out, 
                    const DNAVector & a, int one, 
                    const DNAVector & b, int two) 
{
    out.resize(m_kk);
    int flank = (m_kk - m_k)/2;
        
    int startA = one-flank;
    int stopA = one+m_k;
    int startB = two+m_k;
    int stopB = startB + flank;
        
    if (DEBUG) 
        cerr << "weldableKmer(" << startA << "," << stopA << "); (" << startB << "," << stopB << ")" << endl;
        

    if (startA < 0 || stopB >= b.isize()) {
        out.resize(0);
        if (DEBUG) 
            cerr << "range out of bounds" << endl;
        return;
    }
        
    int i;
    int j = 0;
    for (i=startA; i<stopA; i++) {
        out[j] = a[i];
        j++;
    }
    for (i=startB; i<stopB; i++) {
        out[j] = b[i];
        j++;
    }
        
    if (DEBUG)
        cerr << "\tweld candidate: " << out.AsString() << endl;

}
    
    
bool Welder::Weldable(const DNAVector & a, int one, const DNAVector & b, int two, int thresh, string& welding_kmer, int& weldable_kmer_read_count) 
{
    int i;
    DNAVector d; // stores the required weldabler kmer of length kk
    WeldableKmer(d, a, one, b, two); // constructs teh weldable kmer, stores in (d)
    if (d.isize() == 0)
        return false;
        
        
    int count = m_pTab->GetCount(d, 0); // see if the weldable kmer exists among the reads and get read count.
    weldable_kmer_read_count = count;
        

    if (DEBUG) 
        cerr << "got weldable candidate: " << d.AsString() << " with count: " << weldable_kmer_read_count << endl;
        

    if (count >= thresh) {

        /* 
            if (SimpleHalves(d)) {
            cerr << "Error, halves of weldmer " << d.AsString() << " were found to fail the simple test....  FIXME" << endl;
            exit(3);
            } 
        */
            
        welding_kmer = d.AsString();

        return true;
            
    }
    else
        return false;
}


int Welder::weldmer_count_in_reads(const DNAVector weldmer) 
{
    // weldmer created outside function and provided as input.
                
    int count = m_pTab->GetCount(weldmer, 0); // see if the weldable kmer exists among the reads and get read count.
        
    return(count);
}
