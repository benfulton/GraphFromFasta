// GraphFromFastaTests.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <DNAVector.h>
#include <gtest/gtest.h>
#include <KmerAlignCore.h>
#include <DNAUtil.h>

int _tmain(int argc, _TCHAR* argv[])
{
	::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

TEST(vecDNAVectorTest, Reads) 
{
    vecDNAVector dna;
	string aString = "C:\\Users\\Ben\\Documents\\GIT\\GraphFromFasta\\Tests\\test.fasta";
	dna.Read(aString, false, false, true, 1000000);
	EXPECT_EQ(5000, dna.size());
}

TEST(vecDNAVectorTest, AddData) 
{
	vecDNAVector dna;
	DNAVector v[5];
	v[0].SetFromBases("AAGCTCT");
	v[1].SetFromBases("GTCTGAA");
	v[2].SetFromBases("ATTCGCA");
	v[3].SetFromBases("AAGCTCT");
	v[4].SetFromBases("TCGCACA");

	dna.push_back(v[0]);
	dna.push_back(v[1]);
	dna.push_back(v[2]);
	dna.push_back(v[3]);
	dna.push_back(v[4]);
	KmerAlignCore core;
	core.SetNumTables(2);
	TranslateBasesToNumberExact trans;
	trans.SetSize(5);
	core.SetTranslator(&trans);
	assert( core.Table().size() == 1024);
	core.AddData(dna);
	int b = trans.BasesToNumber(v[0],0);
	int c = trans.BasesToNumber(v[4],0);
	EXPECT_EQ( 2, core.Table()[b].GetData().size());
	EXPECT_EQ( 2, core.Table()[c].GetData().size());
	EXPECT_EQ( 1, core.Table()[183].GetData().size());
}

TEST(DNAUtilTest, compute_entropy) 
{
	DNAVector v;
	v.SetFromBases("AAGCTCT");
	ASSERT_FLOAT_EQ( 1.950212, compute_entropy(v.AsString()));
	v.SetFromBases("A");
	ASSERT_FLOAT_EQ( 0, compute_entropy(v.AsString()));
	v.SetFromBases("AAAA");
	ASSERT_FLOAT_EQ( 0, compute_entropy(v.AsString()));
}
