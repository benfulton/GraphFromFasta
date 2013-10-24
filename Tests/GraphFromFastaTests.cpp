// GraphFromFastaTests.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <DNAVector.h>
#include <gtest/gtest.h>

int _tmain(int argc, _TCHAR* argv[])
{
	::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
  return 0;
}

TEST(FactorialTest, HandlesZeroInput) 
{
    vecDNAVector dna;
	string aString = "C:\\Users\\Ben\\Documents\\GIT\\GraphFromFasta\\Tests\\test.fasta";
	dna.Read(aString, false, false, true, 1000000);
	EXPECT_EQ(5000, dna.size());
}

