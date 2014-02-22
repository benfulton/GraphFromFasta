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

bool is_simple_repeat(const DNAVector& kmer, float MAX_RATIO_INTERNALLY_REPETITIVE) {
    
    // in the future, should just do SW alignment
    // for now, compare all half-kmers
    

    string best_left_kmer;
    string best_right_kmer;
    int best_left_pos = 0;
    int best_right_pos = 0;
    float max_repetitive_ratio = 0;


    int mid_kmer_length = kmer.isize()/2;
    
    for (int i = 0; i < mid_kmer_length; i++) {

        for (int j = i+1; j <= mid_kmer_length; j++) {
        
            int ref_kmer_pos = i;
            int internal_kmer_pos = j;
                    

            // compare half-kmers
            
            int bases_compared = 0;
            int bases_common = 0;
     
            stringstream left_kmer;
            stringstream right_kmer;
            

            while (internal_kmer_pos <= j + mid_kmer_length - 1) {
                bases_compared++;
                if (kmer[ref_kmer_pos] == kmer[internal_kmer_pos]) {
                    bases_common++;
                }
                
                //if (DEBUG) {
                    left_kmer << kmer[ref_kmer_pos];
                    right_kmer << kmer[internal_kmer_pos];
                    //}

                
                ref_kmer_pos++;
                internal_kmer_pos++;
            }
            
            float ratio_same = (float)bases_common/bases_compared;
            if (ratio_same > max_repetitive_ratio) {
                max_repetitive_ratio = ratio_same;
                best_left_kmer = left_kmer.str();
                best_right_kmer = right_kmer.str();
                best_left_pos = i;
                best_right_pos = j;
            }


            if ((! DEBUG) && ratio_same >= MAX_RATIO_INTERNALLY_REPETITIVE) {
                
                // quickest way of exiting this function.  No reason to capture the most repetitive kmer since we're not showing it.
                return(true);
            }
            
        }
           
    }

        
    if (DEBUG) {
        #pragma omp critical
        cout << "# " << kmer.AsString() << " most repetitive kmers:: (" << best_left_pos << "," << best_right_pos << ") " 
             << best_left_kmer << ", " << best_right_kmer << " ratioID: "  << max_repetitive_ratio << endl;
    }

    if (max_repetitive_ratio > MAX_RATIO_INTERNALLY_REPETITIVE) {
        
        if (DEBUG) {
            #pragma omp critical
            cout << "# " << kmer.AsString() << " is internally repetitive: (" << best_left_pos << "," << best_right_pos << ") " 
                 << best_left_kmer << ", " << best_right_kmer << " ratioID: "  << max_repetitive_ratio << endl;
        }
        
        return(true);
    }
    else {
        return(false); // not internally repetitive
    }
}



bool SimpleHalves(const DNAVector & d, float MIN_WELD_ENTROPY, float MAX_RATIO_INTERNALLY_REPETITIVE) {
    int len = d.isize();
    int mid_pos = int(len/2);

    DNAVector left;
    left.resize(mid_pos);
    for (int i = 0; i < mid_pos; i++) {
        left[i] = d[i];
    }
    
    DNAVector right;
    right.resize(len - mid_pos);
    int counter = 0;
    for (int i = mid_pos; i < len; i++) {
        right[counter] = d[i];
        
        counter++;
    }

    if (DEBUG) {
        #pragma omp critical
        {
			cout << "## Left half: " << left.AsString() << ", en: " << compute_entropy(left.AsString());
			cout << "\tRight half: " << right.AsString() << ", en: " << compute_entropy(right.AsString()) << endl;
        }
    }
    
	return (compute_entropy(left.AsString()) < MIN_WELD_ENTROPY 
            || 
			compute_entropy(right.AsString()) < MIN_WELD_ENTROPY
            ||
            is_simple_repeat(left, MAX_RATIO_INTERNALLY_REPETITIVE)
            ||
            is_simple_repeat(right, MAX_RATIO_INTERNALLY_REPETITIVE)

            );
}
        
void Add(vecDNAVector & all, DNAVector & add, int & counter) 
{
    #pragma omp critical
    {
        all.push_back(add);
        counter++;
    }
    
}




