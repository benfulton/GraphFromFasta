#include <string>

struct config
{
	std::string aString;
	bool sStrand;
	std::string readString;
	int num_threads;
	bool REPORT_WELDS;
    double glue_factor;
    int min_glue_required;
	int pooling_kmer_size;
	int welding_kmer_size;
	double min_iso_ratio;
    bool bNoWeld;
    std::string scaffolding_filename;
	float MIN_KMER_ENTROPY;
	float MIN_WELD_ENTROPY;  // min entropy for each half of a welding mer (kk)
	float MAX_RATIO_INTERNALLY_REPETITIVE;

	int MAX_CLUSTER_SIZE;
	int MIN_CONTIG_LENGTH;

	config() : MIN_KMER_ENTROPY(1.3), MIN_WELD_ENTROPY(1.3), MAX_RATIO_INTERNALLY_REPETITIVE(0.85), MAX_CLUSTER_SIZE(100), MIN_CONTIG_LENGTH(24) {}
};

