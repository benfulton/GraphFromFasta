#include "DNAVector.h"
#include "DNAUtil.h"

bool IsSimple(const string & d, float MIN_KMER_ENTROPY) 
{
   return compute_entropy(d) < MIN_KMER_ENTROPY;
}

// used to determine complexity of a kmer
float compute_entropy(const string & kmer) {
    
    map<char,int> char_map;
    
	for (string::const_iterator i = kmer.begin(); i != kmer.end(); ++i) {
        char_map[*i]++;
    }
    
    float entropy = 0;
    
    char nucs[] = { 'G', 'A', 'T', 'C' };
    
    for (char *c = nucs; c != nucs+4; ++c) {
        
        int count = char_map[*c];
        
        float prob = (float)count / kmer.size();
        
        if (prob > 0) {
            float val = prob * log(1/prob)/log((float)2);
            entropy += val;
        }
    }
    
    return(entropy);
}

