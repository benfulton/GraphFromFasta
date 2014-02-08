#ifndef _KMERALIGNCORE_H_
#define _KMERALIGNCORE_H_


#include "DNAVector.h"

//const int KMER_SIZE = 5;
const int KMER_SIZE = 12;

class TranslateBasesToNumberExact
{
 public:

  int BasesToNumber(const DNAVector & b, int from) const {
    int i;
    int num = 0;
    int shift = 1;
    for (i=from; i<from+KMER_SIZE; i++) {
      char v = NucIndex(b[i]);
      if (v == -1 || v == 4)
	return -1;

      num += (int)v * shift;
      shift = (shift << 2);
    }
    //cout << "Translated sequence into " << num << endl;
    return num;
  }

    int GetSize() const {return KMER_SIZE;}
	void SetSize(int n) { throw std::exception(); }

	int GetRange() const {return KMER_SIZE;}
	int GetBoundValue() const { return 1 << (2*KMER_SIZE); }

};


class KmerAlignCoreRecord
{
 public:
  KmerAlignCoreRecord(int contig = -1, int pos = -1, int score = -1) {
    m_contig = contig;
    m_pos = pos;
    m_score = score;
  }

  int GetContig() const {return m_contig;}
  int GetPosition() const {return m_pos;}
  int GetScore() const {return m_score;}

  bool operator < (const KmerAlignCoreRecord & k) const {
    if (m_contig == k.m_contig)
      return (m_pos < k.m_pos);
    return (m_contig < k.m_contig);
  }

  bool operator <= (const KmerAlignCoreRecord & k) const {
    if (m_contig == k.m_contig)
      return (m_pos <= k.m_pos);
    return (m_contig <= k.m_contig);
  }
  bool operator == (const KmerAlignCoreRecord & k) const {
    if (m_contig == k.m_contig)
      return (m_pos == k.m_pos);
    return false;
  }

  KmerAlignCoreRecord & operator =(const KmerAlignCoreRecord & r) {
    m_contig = r.m_contig;
    m_pos = r.m_pos;
    if (r.m_score < m_score)
      m_score = r.m_score;
    return *this;
  }

 private:
  int m_contig;
  int m_pos;
  int m_score;
};


void UniqueSort(std::vector<KmerAlignCoreRecord>& );


class KmerAlignCore
{
public:
  KmerAlignCore();

  void SetMax12Mer(int i) {
    m_max12 = i;
  }

  void AddData(const vecDNAVector & b);

  void SortAll();

  void SetTranslator(TranslateBasesToNumberExact * p) {
    m_pTrans = p;   
  }

  bool GetMatches(std::vector<KmerAlignCoreRecord> & matches, const DNAVector & b, int start = 0);

  const std::vector<KmerAlignCoreRecord> & GetMatchesDirectly(const DNAVector & b, int start = 0);

  void SetNumTables(int i) {
    m_numTables = i;
  }

  int GetWordSize() const {return m_numTables * KMER_SIZE; }

  void SetLookAhead(int i) {m_lookAhead = i;}

  int GetLookAhead() const {return m_lookAhead;}


  void SetLAMaxFreq(int i) {m_lookAheadMaxFreq = i;}

  int NumberOfKmersWithHash( int hash)
  {
	  return m_table[hash].size();
  }
private:
  void MergeSortFilter(std::vector<KmerAlignCoreRecord> & result,
		       const std::vector<KmerAlignCoreRecord> & one,
		       const std::vector<KmerAlignCoreRecord> & two);

  int m_numTables;
  int m_lookAhead;
  int m_max12;

  int m_lookAheadMaxFreq;

  std::map<int, std::vector<KmerAlignCoreRecord>> m_table;
  TranslateBasesToNumberExact * m_pTrans;

};




#endif //_KMERALIGNCORE_H_



