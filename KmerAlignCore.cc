#ifndef FORCE_DEBUG
//#define NDEBUG
#endif

#include <algorithm>
#include "KmerAlignCore.h"

typedef std::vector<KmerAlignCoreRecord> KmerAlignCoreRecordStore;
typedef std::vector<KmerAlignCoreRecordStore> KmerAlignCoreRecordStoreTable;


bool IsRepeat(const NumVector & t, int i, int size, int min)
{
    if (t.size() == 0)
        return false;
    
    int r = 0; 
    int n = 0; 
    for (int j=i; j<i+size; j+=4) {
        if (t[j] >= min)
            r++;
        n++;
    }
    if (r >= n / 2) {
        //cout << "Repeat, skipping!" << endl;
        return true;
    } else {
        return false;
    }
}

void UniqueSort(vector<KmerAlignCoreRecord>& m_data)
{
  sort(m_data.begin(), m_data.end());

  long long i;
  long long k = 0;
  for (i=0; i<m_data.size(); i++) {
    m_data[k] = m_data[i];
    while (i+1<m_data.size() &&  !(m_data[k] < m_data[i+1]))
      i++;
  
    k++;    
  }
  m_data.resize(k);
}

KmerAlignCore::KmerAlignCore() 
{
    m_numTables = 2;    
    m_pTrans = NULL;
    m_lookAhead = 0;
    m_lookAheadMaxFreq = 50000;
    
    m_max12 = 0x7FFFFFFF;
}

void KmerAlignCore::AddData(const vecDNAVector & reads)
{
    cerr << "assigning k-mers..." << endl;

    for (int j=0; j<reads.size(); j++) {
        const DNAVector & read = reads[j];
        if (j % 1000 == 0 || j == reads.size()-1)
            cerr << "\rKmerAlignCore- Contigs: " << j << "   ";
        
		// for each kmer in the read 
        int k=0;
        while (k <= read.size()-KMER_SIZE) {
            int hash_value = m_pTrans->BasesToNumber(read, k);
            if (hash_value >= 0) {
				// j,k are the read index and the kmer index
                m_table[hash_value].push_back(KmerAlignCoreRecord(j,k));
            }
            k++;
        }
    }
    
    cerr << endl << "done!" << endl;
    
}

void vec_sort( pair<const int, vector<KmerAlignCoreRecord>>& s )
{
	sort(s.second.begin(), s.second.end());
}

void KmerAlignCore::SortAll()
{
	for_each( m_table.begin(), m_table.end(), vec_sort);
}

bool KmerAlignCore::GetMatches(vector<KmerAlignCoreRecord> & matches, const DNAVector & b, int start)
{
	int size = m_pTrans->GetSize();
	int i, j, k, l;
	k = start;

	if (start + m_numTables * size > (int)b.size()) {
		cerr << "Error: sequence length=" << b.size() << " and k-kmer end is " << start + m_numTables * size << endl; 
		return false;
	} 


	//cout << "Start position: " << start << endl;

	KmerAlignCoreRecordStoreTable hits;

	hits.resize(m_numTables);

	for (i=0; i<m_numTables; i++) {
		KmerAlignCoreRecordStore & r = hits[i];

		if (m_lookAhead == 0) {
			int n = m_pTrans->BasesToNumber(b, k + i * size);

			if (n >= 0) {
				KmerAlignCoreRecordStore & s = m_table[n];   
				if (s.size() > m_max12) {
					matches.clear();
					return false;
				}

				// Pre-size it!
				r.reserve(s.size());
				int count = 0;
				//cout << "size=" << s.GetNumRecords() << endl;

				for (j=0; j<s.size(); j++) {
					//cout << "Adding single match, pos=" << s.GetRecord(j).GetPosition() << endl;
					KmerAlignCoreRecord const& record = s[j];
					r.push_back(KmerAlignCoreRecord(record.GetContig(), record.GetPosition() - i * size, count));
					//cout << "Done" << endl;
					count++;
					//r.Add(s.GetRecord(j).GetContig(), s.GetRecord(j).GetPosition() - i * size);
				}
			}
		} else {

			int penalty = 0;
			for (l=0; l<m_lookAhead; l++) {
				if (k + i * size + l + size > (int)b.size())
					break;
				//cout << "n=" << k + i * size + l << endl;
				int n = m_pTrans->BasesToNumber(b, k + i * size + l);

				if (n > 0) {
					KmerAlignCoreRecordStore & s = m_table[n];   

					int count = r.size();
					if (l > 0)
						penalty = 1;

					if (l > 0 && count > m_lookAheadMaxFreq) {
						//cout << "***** Count: " << count << endl;
						continue;
					}
					r.resize(count + s.size());
					for (j=0; j<s.size(); j++) {
						//cout << "Adding single match, pos=" << s.GetRecord(j).GetPosition() << endl;
						KmerAlignCoreRecord const& record = s[j];
						r[count] = KmerAlignCoreRecord(record.GetContig(), record.GetPosition() - i * size - l, penalty);
						count++;
						//r.Add(s.GetRecord(j).GetContig(), s.GetRecord(j).GetPosition() - i * size - l);
					}
				}
				//cout << "Before sort: " << r.GetSize() << endl;
				UniqueSort(r);      
				//cout << "After sort:  " << r.GetSize() << endl;
				//for (l=0; l<r.GetSize(); l++) {
				//const KmerAlignCoreRecord & record = r.Get(l);
				//cout << "   c" << record.GetContig() << " @" << record.GetPosition() << endl;
				//}
			}
		}
	}

	matches.clear();
	if (hits.empty())
		return false;

	//cout << "Merging, hits size=" << hits.GetSize() << endl;
	vector<KmerAlignCoreRecord>& tmp = hits[0];
	for (i=1; i<m_numTables; i++) {
		MergeSortFilter(matches, tmp, hits[i]);
		tmp = matches;
		//cout << "Matches remaining: " << matches.isize() << endl;
	}

	if (m_numTables == 1)
		matches = tmp;
	//cout << "All done!" << endl;

	return (matches.size() > 0);
}

void KmerAlignCore::MergeSortFilter(vector<KmerAlignCoreRecord> & result,
	const vector<KmerAlignCoreRecord> & one,
	const vector<KmerAlignCoreRecord> & two) {
		// Stupid for now...
		int i;

		result.clear();

		//cout << "one=" << one.isize() << "  two=" << two.isize() << endl;
		if (one.empty() || two.empty())
			return;

		vector<KmerAlignCoreRecord>::size_type total_size = one.size() + two.size();
		vector<KmerAlignCoreRecord> tmp(total_size);
		result.resize(total_size);

		int k = 0;
		int x = 0;
		int y = 0;
		for (x=0; x<one.size(); x++) {
			while (y<two.size() && two[y] < one[x]) {
				y++;
			}
			if (x >= one.size() || y >= two.size())
				break;
			if (one[x] == two[y]) {
				result[k] = one[x];
				k++;
			}
		}


		result.resize(k);
}

